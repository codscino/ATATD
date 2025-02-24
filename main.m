%% start
tic % the script takes 8 minutes to run on my Apple M1 pro 8 cores
% MATLAB commands
clc, clear, close all, dbstop if error

% Constants and Probe Parameters
muSun   = getAstroConstants('Sun','mu');
muEarth = getAstroConstants('Earth','mu');
RE      = getAstroConstants('Earth', 'Radius');
RS      = getAstroConstants('Sun', 'Radius');
g0      = getAstroConstants('g0');
AU = getAstroConstants('AU');
v_earth = sqrt(muSun/AU); % circular heliocentric earth velocity

% Asteroid Ephemeris Data Loading
data = load('allNEOEphemeris_ATATD_2018_2019.mat');  
neo  = data.allNEOEphemeris_ATATD_2018_2019;

% Append a 9th column with the row index
neo_with_index = [neo, (1:size(neo,1))'];
% Columns:
% 1: semimajor axis [AU]
% 2: eccentricity
% 3: inclination [deg]
% 4: asc. node/raan [deg]
% 5: arg. perigee [deg]
% 6: mean anomaly [deg] at time given in col. 8
% 7: absolute magnitude
% 8: time at which M is given [MJD]
% 9: row index

%% Launch window
% Define departure time window (in MJD2000)
date1_min = date2mjd2000([2033, 01, 01, 0, 0, 0]);
date4_max = date2mjd2000([2037, 12, 31, 0, 0, 0]);

% Define the launch grid (departure dates)
step1 = 15;
% step1 = 5;
grid1 = date1_min:step1:date4_max-500; 
L1 = length(grid1);

% Define the time-of-flight (ToF) grid [days]
ToF1 = 150:step1:450;
L2 = length(ToF1);

%% Asteroid pruning
% FoM filter
e   = neo_with_index(:, 2);               % eccentricity (column 2)
i_deg = neo_with_index(:, 3);           % inclination in degrees (column 3)
i_rad = deg2rad(i_deg);                    % convert inclination to radians

param = sqrt( e.^2 + (2*sin(i_rad/2)).^2 );

% Sort the asteroids based on the computed parameter
[~, sort_idx] = sort(param);

% Extract the top 100 asteroids with the lowest parameter values
top100_indices = sort_idx(1:60);
top100 = neo_with_index(top100_indices, :);

Ltop = size(top100,1);


%% I Batch Conditions
dvData1_cell = cell(Ltop,1);

parfor i = 1:Ltop
    astID = top100(i, 9);
    dvData01 = Leg1ConditionsPar(astID, grid1, ToF1, v_earth, muSun);
    dvData1_cell{i} = dvData01;
end

% Filter asteroids that produced at least one valid transfer.
validFlags = cellfun(@(x) ~isempty(x), dvData1_cell);
idx_validAst = find(validFlags);

% Filter the valid asteroids from top100.
okAsteroids1 = top100(idx_validAst, :);

% Also filter dvData for these asteroids.
dvDataLeg1 = dvData1_cell(idx_validAst);

fprintf(['Number of asteroids passes Leg1' ...
    ' conditions: %d\n'], numel(idx_validAst));

%% II Batch Conditions (Reentry Leg)
numOk = size(okAsteroids1, 1);
dvData2_cell = cell(numOk, 1);

% Define additional parameters for the reentry leg.
step3 = 15;
ToF2 = 50:step3:450;

parfor kf = 1:numOk
    % For each valid asteroid, pass its row from okAsteroids and its
    % corresponding first-stage dv data from dvData1.
    dvData02 = Leg2ConditionsPar(dvDataLeg1{kf}, ...
                  date4_max, step3, ToF2, RE, muEarth, muSun);
    dvData2_cell{kf} = dvData02;
end

% Filter out asteroids that produced no valid reentry entries.
nonEmptyIdx = cellfun(@(x) ~isempty(x), dvData2_cell);
dvData2_cell = dvData2_cell(nonEmptyIdx);
okAsteroids2 = okAsteroids1(nonEmptyIdx, :); % keep only asteroids with valid reentry entries

% Concatenate the results into one final structure array.
validTransfers = vertcat(dvData2_cell{:});
fprintf('Number of asteroids passes both Leg1 and Leg2 conditions: %d \n', ...
    size(okAsteroids2,1));
fprintf('Number of valid transfers: %d\n', numel(validTransfers));


%% final analysis
% Define the finer grid step
stepF = 5;

% Preallocate a cell array for refined results
refinedData = cell(numel(dvData2_cell), 1);

parfor idx = 1:numel(dvData2_cell)
    
    currentData = dvData2_cell{idx};

    % Extract date1 and date2 from the current structure array.
    date1_vals = [currentData.date1_MJD];
    date2_vals = [currentData.date2_MJD];
    
    % Find the minimum date1 in the current cell.
    min_date1 = min(date1_vals);
    max_date1 = max(date1_vals);
    
    % Compute ToF1 for each entry 
    tof1_vals = date2_vals - date1_vals;
    
    % Find the maximum and minimum ToF1 in the current cell.
    max_ToF1 = max(tof1_vals);
    min_ToF1 = min(tof1_vals);
    
    % Calculate the unadjusted starting value:
    rawStart = min_date1 - stepF;
    % Compute the offset from the reference date:
    offset = rawStart - date1_min;
    % Round down the offset to the nearest integer multiple of stepF:
    n = floor(offset / stepF);
    % Recalculate the starting value:
    new_start = date1_min + n * stepF;
    
    % Now define new_grid1 with the adjusted starting value:
    new_grid1 = new_start:stepF:(max_date1 + stepF);

    % Compute new refined grids using stepF 
    % new_grid1 = min_date1-stepF : stepF : max_date1+stepF;
    new_ToF1  = min_ToF1-stepF  : stepF : max_ToF1+stepF;
    
    % Retrieve the asteroid ID from okAsteroids2
    astID = okAsteroids2(idx, 9);
    
    % Re-run Leg1Conditions using the refined grids
    refinedLeg1Data = Leg1ConditionsPar(astID, new_grid1, new_ToF1, v_earth, muSun);

    % Extract date1 and date2 from the current structure array.
    date3_vals = [currentData.date3_MJD];
    date4_vals = [currentData.date4_MJD];

    % Compute ToF1 for each entry 
    tof2_vals = date4_vals - date3_vals;
    
    % Find the maximum and minimum ToF2 in the current cell.
    max_ToF2 = max(tof2_vals);
    min_ToF2 = min(tof2_vals);
    new_ToF2  = min_ToF2-stepF  : stepF : max_ToF2+stepF;

    % Re-run Leg2Conditions using the refined Leg1 output and refined ToF grid
    refinedLeg2Data = Leg2ConditionsPar(refinedLeg1Data, ...
                        date4_max, stepF, new_ToF2, RE, muEarth, muSun);
    
    % Store the refined result for this asteroid.
    refinedData{idx} = refinedLeg2Data;
    
end


%% safe window 15 days reentry
% Copy refinedData into safeData
safeData = refinedData;

for idx = 1:numel(safeData)
    currentData = safeData{idx};
    if isempty(currentData)
        continue;
    end
    
    % Extract date1 and date2 values as column vectors
    date1_vals = [currentData.date1_MJD]';
    date2_vals = [currentData.date2_MJD]';
    
    % Find unique (date1, date2) pairs
    uniquePairs = unique([date1_vals, date2_vals], 'rows');
    
    % Initialize a logical index array to keep rows
    keepRows = true(length(currentData), 1);
    
    % Loop over each unique (date1, date2) pair
    for p = 1:size(uniquePairs, 1)
        u_date1 = uniquePairs(p, 1);
        u_date2 = uniquePairs(p, 2);
        
        % Find indices corresponding to this (date1, date2) pair
        pairIdx = find(date1_vals == u_date1 & date2_vals == u_date2);
        
        % Extract the corresponding date3 values
        date3_vals = [currentData(pairIdx).date3_MJD];
        
        % Compute the spread in date3 values
        date3_spread = max(date3_vals) - min(date3_vals);
        
        % If the spread in date3 is < 15 days, mark these rows for removal
        if date3_spread < 15
            keepRows(pairIdx) = false;
        end
    end
    
    % Update the cell with only the rows that satisfy the condition
    safeData{idx} = currentData(keepRows);
end


%% best dates
% final candiate is astID = 9076, aka 2014 WX202
dvData = safeData{1};

% 1. Min dvTot
dvTot_vec = [dvData.dvTot];
[mindvTot,  minIndexdvTot] = min(dvTot_vec);

% 2. Min dv234
dv234_vec = [dvData.dv234];
[minValue234, minIndex234] = min(dv234_vec);

% 3. Min dv3
dv3_vec = [dvData.dv3];
[minValue3, minIndex3] = min(dv3_vec);

% 4. Min ToF
ToF_vec = [dvData.ToF];
[minValueToF, minIndexToF] = min(ToF_vec);

% 5. Max Layover
lay_vec = [dvData.layover];
[maxValueLay, maxIndexLay] = max(lay_vec);


%%%%%%%%%%%%%%%%%%%%%%%%
%%% 6. Advised date %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Find the minimum date4_MJD over all elements
allDate4 = [dvData.date4_MJD];
minDate4 = min(allDate4);

% Step 2: Select every element where date4_MJD is less than minDate4 + 30
sel1 = dvData(allDate4 < (minDate4 + 30));

% Step 3: Among these, find the maximum date1_MJD
allDate1 = [sel1.date1_MJD];
maxDate1 = max(allDate1);

% Step 4: Further select those elements where maxDate1 - date1_MJD > 15
sel2 = sel1((maxDate1 - [sel1.date1_MJD]) > 15);

% Step 5: Among these, find the element with the lowest dv3
[~, AdvisedIdx] = min([sel2.dv3]);


%% Porkchop plots
selectedIndex = [minIndexdvTot, minIndex234, minIndex3, minIndexToF,...
    maxIndexLay, AdvisedIdx];

porkchopPlotter(dvData, stepF, selectedIndex, 1, 1) % dv1
porkchopPlotter(dvData, stepF, selectedIndex, 1, 2) % dv2
porkchopPlotter(dvData, stepF, selectedIndex, 2, 1) % dv3
porkchopPlotter(dvData, stepF, selectedIndex, 2, 2) % dv4

% Porkchop plots report
% porkchopPlotter_report(dvData, stepF, selectedIndex, 1, 1) % dv1
% porkchopPlotter_report(dvData, stepF, selectedIndex, 1, 2) % dv2
% porkchopPlotter_report(dvData, stepF, selectedIndex, 2, 1) % dv3
% porkchopPlotter_report(dvData, stepF, selectedIndex, 2, 2) % dv4

%% verification results with ODE solver @Advised_date
[posErr, velErr] = odeError(dvData(selectedIndex(6)), muSun);
fprintf('\n'); % Print a blank line
fprintf('Position error ODE: %g km\n', posErr);
fprintf('Velocity error ODE: %g km/s\n', velErr);


%% save results
save('results.mat', 'dvData', 'selectedIndex');
toc
