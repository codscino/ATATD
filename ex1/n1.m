%% start
% matlab commands
tic
clear
clc
close all
dbstop if error

% constants
muSun = getAstroConstants('Sun','mu');
muEarth = getAstroConstants('Earth','mu');
muMars = getAstroConstants('Mars','mu');

RE = getAstroConstants('Earth', 'Radius');
RM = getAstroConstants('Mars', 'Radius');
RS = getAstroConstants('Sun', 'Radius');

% launcher parameter
C3 = [0 5 10 15 20 25 30 40]; % C3 = vinf^2
% C3 = [0 5 6 7 8 9 10 11]; % C3 to make it fail
wetMass = [5750 5000 4250 3600 3000 2400 1900 1000]; % wet mass payload
dryMass = 1250; % dry mass orbiter [kg]
Isp = 310; % Isp orbiter [s]
g0 = getAstroConstants('g0');


%% launch window
%launch window broad range
window1 = date2mjd2000([2026 01 01 0 0 0]);
window2 = date2mjd2000([2027 12 31 12 0 0]);

% find the hohmann transfer date(min theoretical deltav)
[hohmann_date, Thohmann, phasing_target] = phasing_finder(window1,window2);

%define launch window
step1 = 5;
grid1 = hohmann_date-0.67*Thohmann:step1:hohmann_date+0.33*Thohmann;
L1 = length(grid1);

% dumb grid2
step2 = 15;
grid2 = hohmann_date+0.5*Thohmann:step2:hohmann_date+2*Thohmann;
L2 = length(grid2);



%% deltav matrixes for the porkchop
deltavTot = zeros(L1,L2); % array total deltav
deltav1Tot = zeros(L1,L2); % array departure deltav
deltav2Tot = zeros(L1,L2); % array arrival deltav

for i = 1:L1
    for j = 1:L2
        [r1,v1] = EphSS_car(3, grid1(i)); %earth
        [r2,v2] = EphSS_car(4, grid2(j)); %mars WTF
        dT = (grid2(j)-grid1(i)) * 86400;
        
        %%%%%%%%%%%%%%
        %%% tm = 1 %%%
        %%%%%%%%%%%%%%

        [vSc_initial, vSc_final] = myLambert(r1,r2,dT,1,muSun);
        
        % vinf of the spacecraft(earth-centric) 
        v_inf_Earth = norm(vSc_initial-v1); % move to earth-centric
        v_inf_Mars = norm(vSc_final-v2); % move to mars-centric

        % earth departure
        Rc = RE + 250; % parking orbit[km]
        % circular velocity parking orbit (earth centric)
        Vce = sqrt(muEarth/Rc);
        
        % velocity peripasis in escape orbit with vis viva (earth centric)
        Vpe = sqrt(2*(muEarth/Rc+v_inf_Earth^2/2)); 
        
        deltav1_short = Vpe-Vce; % deltaV I manouevre escaping earth
        
        % mars arrival
        Rf = RM+400; %circular arrival orbit around mars radius [km]
        % velocty periapsis arrival hyperbolic orbit to mars [km/s]
        Vpm = sqrt(2*(muMars/Rf+v_inf_Mars^2/2));
        Vcm = sqrt(muMars/Rf); %velocity in operational orbit
        
        deltav2_short = Vpm-Vcm; %deltaV II manouvre to catch Mars
        

        %%%%%%%%%%%%%%%
        %%% tm = -1 %%%
        %%%%%%%%%%%%%%%

        [vSc_initialM, vSc_finalM] = myLambert(r1,r2,dT,-1,muSun);

        % vinf of the spacecraft(earth-centric) 
        v_inf_EarthM = norm(vSc_initialM-v1); % move to earth-centric
        v_inf_MarsM = norm(vSc_finalM-v2); % move to mars-centric
        
        % velocity peripasis in escape orbit with vis viva (earth centric)
        VpeM = sqrt(2*(muEarth/Rc+v_inf_EarthM^2/2)); 
        
        deltav1_long = VpeM-Vce; % deltaV I manouevre escaping earth
        
        % mars arrival

        % velocty periapsis arrival hyperbolic orbit to mars [km/s]
        VpmM = sqrt(2*(muMars/Rf+v_inf_MarsM^2/2)); 
        
        deltav2_long = VpmM-Vcm; %deltaV II manouvre to catch Mars
   
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% final calculations %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        deltavtot_short = deltav1_short + deltav2_short;
        deltavtot_long = deltav1_long + deltav2_long;
        deltavTot(i,j) = min(deltavtot_short,deltavtot_long);
        deltav1Tot(i,j) = min(deltav1_short,deltav1_long);
        deltav2Tot(i,j) = min(deltav2_short,deltav2_long);
    end
end


%% compute maxDV1
clc
% find mindeltav2
minDV2 = min(min(deltav2Tot));

% min wet mass with tsiolkowsky
minWetMass = dryMass*exp(minDV2*1000/(Isp*g0)); %*1000 to convert to m/s
C3max = interp1(wetMass, C3, minWetMass); %C3max si when i have minWetMass
vinfMax = sqrt(C3max);

% velocity peripasis in escape orbit with vis viva (earth centric)
Vpemax = sqrt(2*(muEarth/Rc+vinfMax^2/2));
maxDV1 = norm(Vpemax-Vce);

%% find deltav2min
deltav2Tot_filtered = deltav2Tot;    % Copy original matrix
deltav2Tot_filtered(deltav1Tot > maxDV1) = NaN;  % filtered matrix

foundCandidate = false;  % flag to exit loop once a candidate is accepted

while ~foundCandidate
    % Find the minimum value in deltav2Tot_filtered (ignoring NaNs)
    [currentMin, linearIndex] = min(deltav2Tot_filtered(:));
    
    if isnan(currentMin)
        error(sprintf(['No valid solution found:' ...
            ' All candidates have been rejected.\n' ...
                'You need either a more capable' ...
                ' launcher (greater C3),\n' ...
                'or a more performant orbiter' ...
                ' (higher Isp, lower drymass).']));
    end
    
    % Get row,column indices corresponding to the current minimum candidate
    [i_min, j_min] = ind2sub(size(deltav2Tot_filtered), linearIndex);
    
    % Extract the corresponding departure deltaV candidate from deltav1Tot
    currentDV1 = deltav1Tot(i_min, j_min);
    
    % Invert the Earth departure burn to compute v_inf ---
    vinf_candidate = sqrt( (Vce + currentDV1)^2 - 2*muEarth/Rc );
    
    % Compute C3 from v_inf ---
    C3_candidate = vinf_candidate^2;
    
    % Interpolate to find the corresponding wet mass ---
    wetMass_candidate= interp1(C3, wetMass,C3_candidate,'linear','extrap');
    
    % Compute max arrival deltaV(deltav2_available) with candidate wetMass
    deltav2_available = Isp * g0 * log(wetMass_candidate / dryMass)/1000;
    
    % Check if the available deltaV2 is enough
    if deltav2_available >= currentMin
        % Candidate accepted
        deltav2_min = currentMin;

        % Convert the MJD2000 dates to full date arrays
        start_date = mjd20002date(grid1(i_min));  %%[yyyy,mm,d,h,min,sec]
        arrival_date = mjd20002date(grid2(j_min));
    
        % Get the month names
        start_monthName = monthName(start_date(2));
        arrival_monthName = monthName(arrival_date(2));
    
        fprintf('Minimum deltav2 transfer found at index (%d, %d):\n', ...
            i_min, j_min);
        fprintf('  deltav1 = %.4f km/s\n', currentDV1);
        fprintf('  deltav2 (required) = %.4f km/s\n', currentMin);
        fprintf('  Available deltav2(from wet mass) = %.4f km/s\n', ...
            deltav2_available);
        fprintf('  deltavTot = %.4f km/s\n', currentDV1+currentMin);
        
        % Print the start and arrival dates in a friendly format
        fprintf('  Start date: %s %d, year %d\n', start_monthName, ...
            start_date(3), start_date(1));

        fprintf('  Arrival date: %s %d, year %d\n', arrival_monthName, ...
            arrival_date(3), arrival_date(1));

        transfer_time_dv2min = grid2(j_min) - grid1(i_min);
        fprintf(' Total transfer time: %.1f days\n', transfer_time_dv2min);
    
        foundCandidate = true;
    else
        % Candidate is not acceptable: remove it from and loop again
        deltav2Tot_filtered(i_min, j_min) = NaN;
        fprintf(['Rejected candidate at index (%d, %d): ' ...
            'required deltav2 =' ...
' %.4f km/s but only %.4f km/s available.\n'], ...
            i_min, j_min, currentMin, deltav2_available);
    end
end

% At this point, deltav2_min holds the chosen minimum deltav2 value.

%% Filter dates: Keep only dates (i,j) where both conditions are satisfied:
%  1. deltav1Tot < maxDV1, and
%  2. the available deltav2 (given the launcher wet mass) is
%     more than the required deltav2Tot(by the porkchop)

% start by copying the full matrices into new ones:
deltav1Tot_filtered         = deltav1Tot;
deltav2Tot_filtered_full   = deltav2Tot;  % start from the full matrix
deltavTot_filtered          = deltavTot;

% Loop over all indices and apply the filtering:
for i = 1:L1
    for j = 1:L2
        % First, if deltav1>maxDV1 then mark this index as invalid:
        if deltav1Tot(i,j) > maxDV1
            deltav1Tot_filtered(i,j)       = NaN;
            deltav2Tot_filtered_full(i,j) = NaN;
            deltavTot_filtered(i,j)          = NaN;
        else
            % For indices with deltav1<maxDV1, check the wet mass condition
            % Compute the candidate v_inf using the departure burn:
            vinf_candidate = sqrt( (Vce + deltav1Tot(i,j))^2-2*muEarth/Rc);
            C3_candidate   = vinf_candidate^2;
            wetMass_candidate = interp1(C3, wetMass, ...
                C3_candidate, 'linear', 'extrap');
            % Compute available deltav2 (convert to km/s)
            deltav2_available = Isp*g0*log(wetMass_candidate/dryMass)/1000;
            
            % If the available deltav2 < required deltav2:
            if deltav2_available < deltav2Tot(i,j)
                deltav1Tot_filtered(i,j)       = NaN;
                deltav2Tot_filtered_full(i,j) = NaN;
                deltavTot_filtered(i,j)          = NaN;
            end
        end
    end
end


%% final conclusions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% extent of the launch window %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find the indices (row and column) of non-NaN entries in deltavTot_filtered
[nonNanRows, nonNanCols] = find(~isnan(deltavTot_filtered));

% Determine the minimum and maximum row numbers among these indices
minRow = min(nonNanRows);
maxRow = max(nonNanRows);

window_extent = (maxRow - minRow)*step1;

% Convert the earliest and latest launch dates (MJD2000) to Gregorian date 
earliest_date_array = mjd20002date(grid1(minRow));  %[yyyy,mm,d,h,min,sec]
latest_date_array   = mjd20002date(grid1(maxRow));

% Use your monthName function to get the month name.
earliest_month = monthName(earliest_date_array(2));
latest_month   = monthName(latest_date_array(2));

% Create date strings (without hour, minute, second)
fprintf('\n'); % Print a blank line
earliest_launch = sprintf('%s %d, %d', earliest_month, ...
    earliest_date_array(3), earliest_date_array(1));
latest_launch   = sprintf('%s %d, %d', latest_month, ...
    latest_date_array(3), latest_date_array(1));

% Print the results
fprintf('The launch window extent is %.0f days.\n', window_extent);
fprintf('From %s to %s.\n', earliest_launch, latest_launch);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% baseline launch date %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Since all elements in the deltavTot matrices are feasible launches,
% the safest bet is to launch as early as possible.
bl_date = grid1(minRow);

% Find the indices in the earliest row that are feasible:
feasibleCols = find(~isnan(deltav1Tot_filtered(minRow, :)));
if isempty(feasibleCols)
    error('No feasible baseline candidate in the earliest row.');
end

% Now find the j value in grid2 that gives the minimum time-of-flight 
[~, idx] = min(grid2(nonNanCols) - grid1(minRow));
baselineCol = feasibleCols(idx);

% Extract the baseline candidate parameters from the filtered matrices.
baseline_dV1 = deltav1Tot_filtered(minRow, baselineCol);
baseline_dV2 = deltav2Tot_filtered_full(minRow, baselineCol);
baseline_dVTot = deltavTot_filtered(minRow, baselineCol);

% Convert the baseline start and arrival dates from MJD2000 to Gregorian.
baseline_start_date = mjd20002date(bl_date);% start date from grid1(minRow)
baseline_arrival_date = mjd20002date(grid2(baselineCol)); % arrival date 

% Get month names for printing.
baseline_start_month = monthName(baseline_start_date(2));
baseline_arrival_month = monthName(baseline_arrival_date(2));

% Print the baseline launch candidate details.
fprintf('\n');
fprintf('Baseline launch date found at index (%d, %d):\n', ...
    minRow, baselineCol);
fprintf('  deltav1 = %.4f km/s\n', baseline_dV1);
fprintf('  deltav2 = %.4f km/s\n', baseline_dV2);
fprintf('  deltavTot = %.4f km/s\n', baseline_dVTot);
fprintf('  Start date: %s %d, year %d\n', ...
    baseline_start_month, baseline_start_date(3), baseline_start_date(1));
fprintf('  Arrival date: %s %d, year %d\n', baseline_arrival_month, ...
    baseline_arrival_date(3), baseline_arrival_date(1));

% Compute and print the total transfer time (in days)
transfer_time_min = grid2(baselineCol) - grid1(minRow);
fprintf('  Total transfer time: %.1f days\n', transfer_time_min);

%% Extract dates for hohmann, deltavtotmin and deltav2min
[desired_indexes,desired_minDV] = HdesiredValues(grid1, grid2, ...
    hohmann_date, Thohmann, maxDV1, deltavTot, deltav1Tot, deltav2Tot);

% concatenate deltav2min values
desired_indexes = [desired_indexes; i_min, j_min];
desired_minDV = [desired_minDV;...
    deltavTot(i_min, j_min), deltav1Tot(i_min, j_min), deltav2_min];

% concatenate baseline launch date values
desired_indexes = [desired_indexes; minRow, baselineCol];
desired_minDV = [desired_minDV;...
    baseline_dVTot, baseline_dV1,  baseline_dV2];
%% Porkchop plots
% deltavTot
Hplotmypork(deltavTot_filtered, grid1, grid2, desired_indexes, maxDV1)

% deltav1Tot departure
Hplotmypork(deltav1Tot_filtered, grid1, grid2, desired_indexes, maxDV1)

% deltav2Tot arrival
Hplotmypork(deltav2Tot_filtered_full, grid1, grid2, desired_indexes,maxDV1) 

toc

%% Animated transfer orbit between Earth and Mars        

endTime = grid2(baselineCol);
% simTime = endTime - grid1(1);
 
% figure settings
fig = figure;
set(fig, 'Color', 'k'); 
hsim = axes('Units', 'normalized', 'Position', [0 0.1 1 0.8]);
set(hsim, 'Color', 'none');
set(hsim, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w', 'GridColor', 'w');
set(hsim, 'Layer', 'top');
axes(hsim);

hold on;
title('Transfer Orbit Simulation', 'Color', 'w');
grid on;
axis equal;

% set image as background of the plot (comment to have simple black backgr)
ha = axes('Units', 'normalized', 'Position', [0 0 1 1]);
I = imread('deepfield4k.tiff');
hi = imagesc(I);
set(ha, 'HandleVisibility', 'off', 'Visible', 'off');
uistack(ha, 'bottom');

% plot 3d Sun static in the middle
plot_planet(0, 0, 0, RS * 30, 'ATATD-Toolbox/TButils/textures/Sun.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% plot planets orbit %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate periods
T_earth = getAstroConstants('Earth', 'Period') / 86400;
T_mars  = getAstroConstants('Mars', 'Period') / 86400;

% calculate initial poistions of earth and mars
[r1e, v1e] = EphSS_car(3, grid1(1));  
[r1m, v1m] = EphSS_car(4, grid1(1)); 

% plot them for their period + 1 day(to be sure ellipse closes)
nPointsEarth = ceil(T_earth) + 1;
nPointsMars  = ceil(T_mars) + 1;
r_earth_static = zeros(nPointsEarth, 3);
r_mars_static  = zeros(nPointsMars, 3);

for i = 1:nPointsEarth
    r_earth_static(i, :) = FGKepler_dt2(r1e, v1e, i * 86400, muSun);
end
for i = 1:nPointsMars
    r_mars_static(i, :) = FGKepler_dt2(r1m, v1m, i * 86400, muSun);
end

plot3(r_earth_static(:,1), r_earth_static(:,2), ...
    r_earth_static(:,3), 'b', 'LineWidth', 2);
plot3(r_mars_static(:,1), r_mars_static(:,2), ...
    r_mars_static(:,3), 'r', 'LineWidth', 2);


% initialize 3d planets models
hEarthModel = plot_planet(r1e(1), r1e(2), r1e(3), RE * 3000, ...
    'ATATD-Toolbox/TButils/textures/Earth.jpg');
hMarsModel  = plot_planet(r1m(1), r1m(2), r1m(3), RM * 3000, ...
    'ATATD-Toolbox/TButils/textures/Mars.jpg');

% initialize counters text labels
hDateText = text(0.05, 0.95, '', 'Units', 'normalized', ...
    'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w');
hElapsedText = text(0.05, 0.88, '', 'Units', 'normalized', ...
    'FontSize', 12, 'Color', 'm', 'BackgroundColor', 'w');
hPhasingText = text(0.05, 0.82, '', 'Units', 'normalized', ...
    'FontSize', 12, 'Color', 'k', 'BackgroundColor', 'w');

% initialize phasing angle lines
hLineEarth = plot3(nan, nan, nan, 'w-', 'LineWidth', 0.5);
hLineMars  = plot3(nan, nan, nan, 'w-', 'LineWidth', 0.5);
hPhasingArc = plot3(nan, nan, nan, 'w-', 'LineWidth', 0.5);

% initialize slider
hSlider = uicontrol('Style', 'slider', 'Min', grid1(1), ...
    'Max', endTime, 'Value', grid1(1), ...
    'Units', 'normalized', 'Position', [0.2, 0.02, 0.6, 0.05]);

% actual transfer orbit calculations
t_depart = grid1(minRow);
t_arrive = grid2(baselineCol);
dT_transfer = (t_arrive - t_depart) * 86400;
[r_depart, v_depart] = EphSS_car(3, t_depart);
[r_arrive, v_arrive] = EphSS_car(4, t_arrive);
% baseline launch date tm=-1(long arc) because period longer than Thohmann
[vSc_initial_bl, vSc_final_bl] = myLambert(r_depart, r_arrive, ...
    dT_transfer, -1, muSun);
n_points_sc = 1000;
r_sc_traj = zeros(n_points_sc, 3);
for k = 1:n_points_sc
    t_offset = (k - 1) / (n_points_sc - 1) * dT_transfer;
    r_sc_traj(k, :) = FGKepler_dt2(r_depart, vSc_initial_bl, ...
        t_offset, muSun);
end

% initialize transfer orbit and spacecraft position plot
hTransferOrbit = plot3(nan, nan, nan, 'm--', 'LineWidth', 2);
hSpacecraft = plot3(nan, nan, nan, 'o', 'MarkerSize', 8, ...
    'MarkerFaceColor', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slider functions to update the plot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% listener
addlistener(hSlider, 'Value', 'PostSet', @(src, event) ...
    updatePlot(get(hSlider, 'Value'), hDateText, hElapsedText, ...
    hPhasingText, hPhasingArc, hLineEarth, hLineMars, ...
    hEarthModel, hMarsModel, hTransferOrbit, hSpacecraft, ...
    r1e, v1e, r1m, v1m, r_depart, vSc_initial_bl, r_sc_traj, ...
    t_depart, t_arrive, muSun, grid1, RE, RM, phasing_target));

% updatePlot
updatePlot(get(hSlider, 'Value'), hDateText, hElapsedText, ...
    hPhasingText, hPhasingArc, hLineEarth, hLineMars, ...
    hEarthModel, hMarsModel, hTransferOrbit, hSpacecraft, ...
    r1e, v1e, r1m, v1m, r_depart, vSc_initial_bl, r_sc_traj, ...
    t_depart, t_arrive, muSun, grid1, RE, RM, phasing_target);

% updatePlot function
function updatePlot(currentMJD, hDateText, hElapsedText, ...
    hPhasingText, hPhasingArc, hLineEarth, hLineMars, ...
    hEarthModel, hMarsModel, hTransferOrbit, hSpacecraft, ...
    r1e, v1e, r1m, v1m, r_depart, vSc_initial, r_sc_traj, ...
    t_depart, t_arrive, muSun, grid1, ~, ~, hohmannPhasingTarget)
    
    drawnow limitrate;
    
    dateVec = mjd20002date(currentMJD);
    dateStr = sprintf('%d %s %d', dateVec(3), monthName(dateVec(2)), ...
        dateVec(1));
    set(hDateText, 'String', dateStr);
    if currentMJD >= t_depart
        daysElapsed = currentMJD - t_depart;
        set(hElapsedText, 'String', ...
            sprintf('Days since transfer start: %.1f', daysElapsed));
    else
        set(hElapsedText, 'String', '');
    end
    
    offsetPlanets = (currentMJD - grid1(1)) * 86400;
    newEarthCenter = FGKepler_dt2(r1e, v1e, offsetPlanets, muSun);
    newMarsCenter  = FGKepler_dt2(r1m, v1m, offsetPlanets, muSun);
    
    baseEarth = get(hEarthModel, 'UserData');
    newX_Earth = baseEarth.baseX + newEarthCenter(1);
    newY_Earth = baseEarth.baseY + newEarthCenter(2);
    newZ_Earth = baseEarth.baseZ + newEarthCenter(3);
    set(hEarthModel, 'XData', newX_Earth, 'YData', newY_Earth, ...
        'ZData', newZ_Earth);
    
    baseMars = get(hMarsModel, 'UserData');
    newX_Mars = baseMars.baseX + newMarsCenter(1);
    newY_Mars = baseMars.baseY + newMarsCenter(2);
    newZ_Mars = baseMars.baseZ + newMarsCenter(3);
    set(hMarsModel, 'XData', newX_Mars, 'YData', newY_Mars, ...
        'ZData', newZ_Mars);
    
    set(hLineEarth, 'XData', [0, newEarthCenter(1)], 'YData', ...
        [0, newEarthCenter(2)], 'ZData', [0, newEarthCenter(3)]);
    set(hLineMars,  'XData', [0, newMarsCenter(1)],  'YData', ...
        [0, newMarsCenter(2)],  'ZData', [0, newMarsCenter(3)]);
    
    lambdaEarth = atan2d(newEarthCenter(2), newEarthCenter(1));
    lambdaMars  = atan2d(newMarsCenter(2), newMarsCenter(1));
    phaseAngle  = mod(lambdaMars - lambdaEarth, 360);
    
    dangle = abs(phaseAngle - hohmannPhasingTarget);
    if dangle > 180
        dangle = 360 - dangle;
    end
    
    % highlight in green when the angle is +-5 degree of hohmann ph angle
    if dangle <= 5
        set(hPhasingText, 'Color', 'g');
        set(hPhasingArc, 'Color', 'g', 'LineWidth', 1.5);
    else
        set(hPhasingText, 'Color', 'k');
        set(hPhasingArc, 'Color', 'w', 'LineWidth', 0.5);
    end
    set(hPhasingText, 'String', sprintf('Phasing angle: %.1f°', ...
        phaseAngle));
    
    arcRadius = 0.4 * min(norm(newEarthCenter), norm(newMarsCenter));
    startAngle = deg2rad(lambdaEarth);
    endAngle   = startAngle + deg2rad(phaseAngle);
    thetaArc = linspace(startAngle, endAngle, 50);
    xArc = arcRadius * cos(thetaArc);
    yArc = arcRadius * sin(thetaArc);
    zArc = zeros(size(thetaArc));
    set(hPhasingArc, 'XData', xArc, 'YData', yArc, 'ZData', zArc);
    
    if currentMJD < t_depart
        set(hTransferOrbit, 'XData', nan, 'YData', nan, 'ZData', nan);
        set(hSpacecraft, 'XData', nan, 'YData', nan, 'ZData', nan);
    elseif currentMJD >= t_depart && currentMJD <= t_arrive
        offset = (currentMJD - t_depart) * 86400;
        r_sc_current = FGKepler_dt2(r_depart, vSc_initial, offset, muSun);
        N_line = 20;
        t_line = linspace(0, offset, N_line);
        r_line = zeros(N_line, 3);
        for j = 1:N_line
            r_line(j, :) = FGKepler_dt2(r_depart, vSc_initial, ...
                t_line(j), muSun);
        end
        set(hTransferOrbit, 'XData', r_line(:,1), ...
            'YData', r_line(:,2), 'ZData', r_line(:,3));
        set(hSpacecraft, 'XData', r_sc_current(1), ...
            'YData', r_sc_current(2), 'ZData', r_sc_current(3));
    else
        set(hTransferOrbit, 'XData', r_sc_traj(:,1), ...
            'YData', r_sc_traj(:,2), 'ZData', r_sc_traj(:,3));
        set(hSpacecraft, 'XData', r_arrive(1), ...
            'YData', r_arrive(2), 'ZData', r_arrive(3));
    end
    
    drawnow;
end

%% verification with ODE
F = @(t,x) [x(4);               % dx/dt = Vx
             x(5);              % dy/dt = Vy
             x(6);              % dz/dt = Vz
            -muSun*x(1)/((x(1)^2+x(2)^2+x(3)^2)^(3/2));   % dVx/dt
            -muSun*x(2)/((x(1)^2+x(2)^2+x(3)^2)^(3/2));   % dVy/dt
            -muSun*x(3)/((x(1)^2+x(2)^2+x(3)^2)^(3/2) )]; % dVz/dt

% seconds from departure to arrival
tspan = [0, (t_arrive-t_depart)*86400]; 
% initial state (position and velocity from Lambert solver)
y0 = [r_depart, vSc_initial_bl];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
[t,y] = ode45(F, tspan, y0, options);

fpos = y(end, 1:3); % final postion
vpos = y(end, 4:6); % final velocity

%%%%%%%%%%%%%%%%%%
%%% ODE errors %%%
%%%%%%%%%%%%%%%%%%

% position
posErr = norm(r_arrive-fpos);
% rSc_final = FGKepler_dt2(r_depart, vSc_initial_bl,  dT_transfer, muSun);
% posErr = norm(rSc_final-fpos);

% velocity
% account deltav2 braking in the direction of v_arrive
baseline_dV2_vec = -baseline_dV2 * (v_arrive / norm(v_arrive));
% Compute the corrected target arrival velocity (after the braking burn)
v_target = v_arrive + baseline_dV2_vec;
% Now compute the velocity error relative to the target velocity
velErr1 = norm(v_target-vpos);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% velErr2: uses lambert arc final velocity instead ODE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vrel = vSc_final_bl-v_arrive;
velErr2 = abs(norm(vrel)- Vcm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print results with appropriate units %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n'); % Print a blank line
fprintf('Position error ODE (km): %g\n', posErr);
fprintf(['Velocity error 1 (km/s)' ...
    ' [ODE corrected for braking]: %g\n'], velErr1);
fprintf('Velocity error 2 (km/s) [Lambert arc vs. Vcm]: %g\n', velErr2);

%% verify lambert solver using SOTA lambert solver
% Time the SOTA Lambert solver call
tic
[V1, V2, extremal_distances, exitflag] = izzoLambert(r_depart, r_arrive, ...
    -dT_transfer/86400, 0, muSun);
timeLambert = toc;

% Time the gigaFun Lambert solver call
tic
[vSc_initial_bl, vSc_final_bl] = myLambert(r_depart, r_arrive, ...
    dT_transfer, -1, muSun);
timeGigaFun = toc;

% Compute absolute velocity errors (in km/s)
vInErr = norm(V1 - vSc_initial_bl);
vFinErr = norm(V2 - vSc_final_bl);

% Compute percentage difference in computation time relative to gigaFun:
timeDiffPct = (abs(timeLambert - timeGigaFun) / timeGigaFun) * 100;

% Print results
fprintf('\n'); % Print a blank line
fprintf(['My lambert solver ' ...
    'compared to state-of-the-art Lambert solver:\n']);
fprintf('Initial velocity error: %g km/s\n', vInErr);
fprintf('Final velocity error:   %g km/s\n', vFinErr);
fprintf(['Percentage difference ' ...
    'in computation time: %.2f%%\n\n'], timeDiffPct);
