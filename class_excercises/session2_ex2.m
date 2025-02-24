%% start
clear
clc
close all

dbstop if error

muSun = getAstroConstants('Sun','mu');

% dates
date1 = date2mjd2000([2016 03 14 0 0 0]);
date2 = date2mjd2000([2016 10 15 0 0 0]);

% Generate grid of conditions
Range = 3*30; %90 days
step = 5; % 5 days;
grid1 = (date1-Range):step:(date1+Range);
grid2 = (date2-Range):step:(date2+Range);
L = size(grid1,2);

%% planet position grid
r1 = zeros(1,3); % initialize earth positions array
r2 = zeros(1,3); % initialize mars positions array
v1 = zeros(1,3); % initialize earth velocity array
v2 = zeros(1,3); % initialize mars velocity array

vSc_initial = zeros(1,3); % initialize  
vSc_final = zeros(1,3); % initialize 

vSc_initialM = zeros(1,3); % initialize  
vSc_finalM = zeros(1,3); % initialize 

deltavTot = zeros(L,L); % initialize 

for i = 1:L
    for j = 1:L
        [r1,v1] = EphSS_car(3, grid1(i)); %earth
        [r2,v2] = EphSS_car(4, grid2(j)); %mars WTF
        dT = (grid2(j)-grid1(i)) * 86400;

        % tm = 1
        [vSc_initial, vSc_final] = myLambert(r1,r2,dT,1,muSun);
        deltav1 = norm(v1-vSc_initial);
        deltav2 = norm(v2-vSc_final);
        
        % tm = -1
        [vSc_initialM, vSc_finalM] = myLambert(r1,r2,dT,-1,muSun);
        deltav1M = norm(v1-vSc_initialM);
        deltav2M = norm(v2-vSc_finalM);

        deltavtot = deltav1 + deltav2;
        deltavtotM = deltav1M + deltav2M;
        deltavtotMin = min(deltavtot,deltavtotM);
        deltavTot(i,j) = deltavtotMin;
    end
end


%% plot
minDV = min(min(deltavTot));
figure
hold on

% porkchop plot
contourf(grid1, grid2, deltavTot', minDV:0.15:10)
% why do I have to make the transpose here?

% plot minDV point
plot(date1, date2, 'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0], ...
    'MarkerSize', 10, ...
    'Marker', 'square', ...
    'LineStyle','none');


xlabel('MJD departure');
ylabel('MJD arrival');
title('Pork-Chop plot total dV');

% color map from nasa report
colors = [
    255,0,255;
    255,0,0;
    255,102,0;
    255,153,0;
    255,204,0;
    255,255,0;
    0,255,0;
    153,204,0;
    51,204,204;
    0,255,255;
    0,204,255;
    51,102,255;
    0,0,255;
    153,52,102;
    128,0,128]/255; % normalizibg rgb

colormap(colors)
colorbar

% grid part
axis equal
xGrid = linspace(min(grid1(:,:)), max(grid1(:,:)), 20);
yGrid = linspace(min(grid2(:,:)), max(grid2(:,:)), 20);


% Add vertical lines
for x = xGrid
    xline(x, 'k-', 'LineWidth', 1);
end

% Add horizontal lines
for y = yGrid
    yline(y, 'k-', 'LineWidth', 1);
end


% plot same dT lines

%Define delta time values
% delta_t_values = [50, 100, 150, 200, 250, 300, 350, 400];
% 
% for dt = delta_t_values
%     arrival_dates = grid1 + dt; % Compute arrival dates
%     plot(grid1, arrival_dates, 'w', 'LineWidth', 2); % Plot white line
% 
%     % Position label at the middle of each line, slightly to the right
%     mid_index = floor(length(grid1) / 2);
%     x_label = grid1(mid_index) + 5; % Shift slightly to the right
%     y_label = arrival_dates(mid_index);
% 
%     text(x_label, y_label, sprintf('%d days', dt), ...
%         'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', ...
%         'BackgroundColor', 'k', ...
%         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
% end


% Predefined axis limits
x_limits = [5826.5, 6006.5];
y_limits = [6041.5, 6221.5];
axis([x_limits, y_limits]);
