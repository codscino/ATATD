function Hplotmypork(deltavTot, grid1, grid2, desired_indexes, maxDV1, color1, color2, color3, color4, colors)
    
    if nargin < 6
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
        
        % Default colors for the markers:
        % Hohmann date
        color1 = [0, 0, 0];
        % Optimal date
        color2 = [255, 255, 255] / 255;
        % Optimal dV2 (dv1max < maxDV1)
        color3 = [0, 0, 255] / 255;
        % Baseline launch date
        color4 = [0, 255, 153] / 255;
    end

    % Find minimum dV value for contour plotting
    minDV = min(min(deltavTot));
    
    figure;
    hold on;
    axis equal;

    % Define custom tick positions 
    xGrid = linspace(min(grid1(:)), max(grid1(:)), 20);
    yGrid = linspace(min(grid2(:)), max(grid2(:)), 20);
    
    % Plot the porkchop contour
    contourf(grid1, grid2, deltavTot', linspace(minDV, minDV * 1.5, 10));
    
    % =======================================================================
    % Extract the (departure, arrival) points for the four markers.
    % (Assuming desired_indexes is a 4-by-2 matrix with indices into grid1 and grid2)
    % =======================================================================
    points = zeros(4, 2);  % each row: [x, y]
    for k = 1:4
        points(k, :) = [grid1(desired_indexes(k, 1)), grid2(desired_indexes(k, 2))];
    end

    % =======================================================================
    % Plot each marker with different shape, size, and transparency
    % =======================================================================
    h1 = scatter(points(1, 1), points(1, 2), 50, 'o', 'filled', ...
        'MarkerFaceColor', color1, 'MarkerEdgeColor', [0, 0, 0]);
    h1.MarkerFaceAlpha = 1;
    
    h2 = scatter(points(2, 1), points(2, 2), 50, 'o', 'filled', ...
        'MarkerFaceColor', color2, 'MarkerEdgeColor', [0, 0, 0]);
    h2.MarkerFaceAlpha = 1;
    
    h3 = scatter(points(3, 1), points(3, 2), 50, 'o', 'filled', ...
        'MarkerFaceColor', color3, 'MarkerEdgeColor', [0, 0, 0]);
    h3.MarkerFaceAlpha = 1;

    h4 = scatter(points(4, 1), points(4, 2), 50, 'o', 'filled', ...
        'MarkerFaceColor', color4, 'MarkerEdgeColor', [0, 0, 0]);
    h4.MarkerFaceAlpha = 1;
    
    % =======================================================================
    % Labels, title, colormap, and colorbar
    % =======================================================================
    xlabel('Date departure');
    ylabel('Date arrival');
    
    varName = inputname(1);
    if contains(varName, '1')
        plotTitle = 'DeltaV 1';
    elseif contains(varName, '2')
        plotTitle = 'DeltaV 2';
    elseif contains(lower(varName), 'delta')
        plotTitle = 'DeltaV total';
    else
        error(['Invalid inputname variable: "', varName, '". ', ...
               'The variable name must contain either ''1'', ''2'', or ''delta''.']);
    end
    title(['PorkChop plot: ', plotTitle]);
    
    colormap(colors);
    colorbar;
    axis equal;
    
    legend([h1, h2, h3, h4], ...
        {'Hohmann Date', 'Optimal Date (deltavTot min)', ...
        ['Optimal dV2 (dv1max < ', num2str(maxDV1), ')'], ...
        'Baseline launch date (earliest date with minDt)'}, ...
        'Location', 'bestoutside');

    % Set the tick positions and convert to Gregorian date labels
    set(gca, 'XTick', xGrid, 'YTick', yGrid);
    xtickLabels = cell(size(xGrid));
    for k = 1:length(xGrid)
        dateArr = mjd20002date(xGrid(k));
        xtickLabels{k} = sprintf('%d %s', dateArr(3), monthName(dateArr(2)));
    end
    set(gca, 'XTickLabel', xtickLabels);
    
    ytickLabels = cell(size(yGrid));
    for k = 1:length(yGrid)
        dateArr = mjd20002date(yGrid(k));
        ytickLabels{k} = sprintf('%d %s', dateArr(3), monthName(dateArr(2)));
    end
    set(gca, 'YTickLabel', ytickLabels);
    
    % -----------------------------------------------------------------------
    % Instead of manually drawing vertical/horizontal grid lines,
    % we simply turn the built-in grid on and force it to the bottom.
    grid on;
    set(gca, 'Layer', 'bottom');  % Draw grid lines behind other objects.
    % -----------------------------------------------------------------------
    
    % =======================================================================
    % Draw lines of constant transfer time (deltaTime)
    % =======================================================================
    % Transfer time: deltaTime = arrival - departure = y - x.
    % We draw lines y = x + dt for dt = 120, 150, ..., 480 days.
    xLimits = get(gca, 'XLim');
    yLimits = get(gca, 'YLim');
    
    for dt = 120:30:480
        x_start = max(xLimits(1), yLimits(1) - dt);
        x_end   = min(xLimits(2), yLimits(2) - dt);
        if x_start < x_end
            x_line = [x_start, x_end];
            y_line = x_line + dt;
            % Draw the constant dt line (suppress from legend)
            plot(x_line, y_line, 'k-', 'LineWidth', 0.25, 'HandleVisibility', 'off');
            
            % Place a text label at 1/4 of the x-range from x_start.
            x_label = x_start + (x_end - x_start) / 4;
            y_label = x_label + dt;
            text(x_label, y_label, sprintf('%d days', dt), 'Color', 'k', 'FontSize', 7, ...
                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                 'BackgroundColor', 'w');
        end
    end

    hold off;
end