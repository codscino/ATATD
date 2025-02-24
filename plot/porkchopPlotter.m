function porkchopPlotter(dvData, stepF, selectedIndex, legNum, dvNum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% porkchopPlotter - Generate porkchop plots for transfer maneuver data
%
%   Inputs:
%       dvData       - Structure array containing transfer maneuver
%                      entries with date fields and delta-v values.
%       stepF        - Grid refinement step size for plotting (in days).
%       selectedIndex- Array of indices specifying selected entries
%                      for highlighting.
%       legNum       - Integer (1 or 2) specifying which leg of
%                      the transfer to plot.
%       dvNum        - Integer (1 or 2) specifying which delta-v
%                      value to plot for the selected leg.
%
%   Outputs:
%       (No explicit output; the function generates a figure displaying
%        the porkchop plot.)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    

    % Color map from NASA report (normalized RGB values)
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
        128,0,128] / 255;

    if legNum == 1
        date1 = 'date1_MJD';
        date2 = 'date2_MJD';
         if dvNum == 1
            dvnum = 'dv1';
        elseif dvNum == 2
            dvnum = 'dv2';
         end
    elseif legNum == 2
        date1 = 'date3_MJD';
        date2 = 'date4_MJD';
        if dvNum == 1
            dvnum = 'dv3';
        elseif dvNum == 2
            dvnum = 'dv4';
        end
    else
        error('invalid values, please refer to the function documentation')
    end

    minDate1 = min([dvData.(date1)]);
    maxDate1 = max([dvData.(date1)]);
    minDate2 = min([dvData.(date2)]);
    maxDate2 = max([dvData.(date2)]);

    % Define grid vectors using the refined step
    plot_grid1 = minDate1:stepF:maxDate1;
    plot_grid2 = minDate2:stepF:maxDate2;
    rowNames = string(plot_grid1);
    colNames = string(plot_grid2);
    
    % Initialize a table with NaNs (or some other default value)
    dv_table = array2table(nan(numel(rowNames), numel(colNames)), ...
        'RowNames', rowNames, 'VariableNames', colNames);
    
    % Populate the table by assigning values using the date keys
    for j = 1:length(dvData)
        dep_date = dvData(j).(date1);
        arr_date = dvData(j).(date2);
        
        % Convert numeric dates to string keys
        rowKey = string(dep_date);
        colKey = string(arr_date);
        
        % Assign the dv1 value
        dv_table{rowKey, colKey} = dvData(j).(dvnum);
    end
    
    % Convert the table to a numeric matrix.
    dv_grid = table2array(dv_table);
    
    minDV = min(min(dv_grid));
    maxDV = max(max(dv_grid));
    
    
    %Plot the porkchop contour
    figure;
    hold on;
    grid on;
    axis equal;
    % Define contour levels
    contourLevels = linspace(minDV, maxDV, 30);
    contourf(plot_grid1, plot_grid2, dv_grid', ...
        contourLevels, 'LineStyle', 'none');
    colormap(colors);
    cb = colorbar;
    cb.Label.String = 'deltaV [km/s]';
    cb.Label.FontSize = 20;
    xlabel('Departure date','FontSize',24);
    ylabel('Arrival date','FontSize',24);
    title(['Porkchop Plot ' dvnum],'FontSize',20);


    % =====================================================================
    % Plot selected points
    % =====================================================================

    % 1. Min dvTot: Black circle with black edge
    date1_I = dvData(selectedIndex(1)).(date1);
    date2_I = dvData(selectedIndex(1)).(date2);
    h1 = scatter(date1_I, date2_I, 10^2, '^', 'MarkerEdgeColor','k',...
        'MarkerFaceColor','k','MarkerFaceAlpha',1);
    
    % 2. Min dv234: Red square with black edge
    date1_II = dvData(selectedIndex(2)).(date1);
    date2_II = dvData(selectedIndex(2)).(date2);
    h2 = scatter(date1_II, date2_II, 20^2, 's', 'MarkerEdgeColor','k',...
        'MarkerFaceColor','r','MarkerFaceAlpha', 0.8);
    
    % 3. Min dv3: Blue diamond with black edge
    date1_III = dvData(selectedIndex(3)).(date1);
    date2_III = dvData(selectedIndex(3)).(date2);
    h3 = scatter(date1_III, date2_III, 15^2, 'd', 'MarkerEdgeColor','k',...
        'MarkerFaceColor','b','MarkerFaceAlpha', 0.7);
    
    % 4. Min ToF: Green triangle up with black edge
    date1_IV = dvData(selectedIndex(4)).(date1);
    date2_IV = dvData(selectedIndex(4)).(date2);
    h4 = scatter(date1_IV, date2_IV, 10^2, '^', 'MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerFaceAlpha', 0.7);
    
    % 5. Max Layover: Cyan plus sign with thicker line (plus signs are not filled)
    date1_V = dvData(selectedIndex(5)).(date1);
    date2_V = dvData(selectedIndex(5)).(date2);
    h5 = scatter(date1_V, date2_V, 8^2, '+', 'MarkerEdgeColor','c',...
        'LineWidth', 2, 'MarkerFaceAlpha', 0.1);
    
    % 6. Advised date: Magenta pentagram with black edge
    date1_VI = dvData(selectedIndex(6)).(date1);
    date2_VI = dvData(selectedIndex(6)).(date2);
    h6 = scatter(date1_VI, date2_VI, 10^2, 'p', 'MarkerEdgeColor','k',...
        'MarkerFaceColor','m','MarkerFaceAlpha', 0.8);
        
    
    % =====================================================================
    % Set axes labels in gregorian dates
    % =====================================================================
    % Define custom tick positions 
    xGrid = linspace(min(plot_grid1(:)), max(plot_grid1(:)), 5);
    yGrid = linspace(min(plot_grid2(:)), max(plot_grid2(:)), 5);
    % Set the tick positions and convert to Gregorian date labels
    % Set the axes font size for bigger tick labels
    set(gca, 'FontSize', 20);

    set(gca, 'XTick', xGrid, 'YTick', yGrid);
    xtickLabels = cell(size(xGrid));
    for k = 1:length(xGrid)
        dateArr = mjd20002date(xGrid(k));
        xtickLabels{k} = sprintf('%d %s %d', dateArr(3), ...
            monthName(dateArr(2)), dateArr(1));
    end
    set(gca, 'XTickLabel', xtickLabels);

    ytickLabels = cell(size(yGrid));
    for k = 1:length(yGrid)
        dateArr = mjd20002date(yGrid(k));
        ytickLabels{k} = sprintf('%d %s %d', dateArr(3), ...
            monthName(dateArr(2)), dateArr(1));
    end
    set(gca, 'YTickLabel', ytickLabels);




    % =====================================================================
    % Draw lines of constant transfer time (deltaTime)
    % =====================================================================
    xLimits = get(gca, 'XLim');
    yLimits = get(gca, 'YLim');
    
    for dt = 150:30:400
        x_start = max(xLimits(1), yLimits(1) - dt);
        x_end   = min(xLimits(2), yLimits(2) - dt);
        if x_start < x_end
            x_line = [x_start, x_end];
            y_line = x_line + dt;
            % Draw the constant dt line (suppress from legend)
            plot(x_line, y_line, 'k-', 'LineWidth', 0.25, ...
                'HandleVisibility', 'off');
            
            % Place a text label at 1/4 of the x-range from x_start.
            x_label = x_start + (x_end - x_start) / 2;
            y_label = x_label + dt;
            text(x_label, y_label, sprintf('%d days', dt), ...
                'Color', 'k', 'FontSize', 20, ...
                 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', ...
                 'BackgroundColor', 'w');
        end
    end


    % =====================================================================
    % Legend
    % =====================================================================
    legend([h1, h2, h3, h4, h5, h6], ...
        {'Min dvTot','Min dv234','Min dv3', ...
        'Min ToF','Max Layover','Advised date'}, 'Location', 'NorthWest');

   
    
    hold off;
    
end