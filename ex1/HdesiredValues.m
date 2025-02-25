function [desired_indexes,desired_minDV] = HdesiredValues(grid1, grid2, hohmann_date, Thohmann, maxDV1, deltavTot, deltav1Tot, deltav2Tot)
    
    % --- Find the closest hohmann_date index in grid1 ---
    [~, i_h] = min(abs(grid1 - hohmann_date));  % find index where grid1 is closest to hohmann_date
    if abs(grid1(i_h) - hohmann_date) > 10
        warning('Closest hohmann date in grid1 (%g) is more than 10 days away from the target (%g).', grid1(i_h), hohmann_date);
    end
    
    % --- Find the closest (hohmann_date + Thohmann) index in grid2 ---
    target_date = hohmann_date + Thohmann;
    [~, j_h] = min(abs(grid2 - target_date));  % find index where grid2 is closest to (hohmann_date+Thohmann)
    if abs(grid2(j_h) - target_date) > 20
        warning('Closest (hohmann_date + Thohmann) in grid2 (%g) is more than 20 days away from the target (%g).', grid2(j_h), target_date);
    end

    minDV_h = deltavTot(i_h, j_h);
    minDV1_h = deltav1Tot(i_h, j_h);
    minDV2_h = deltav2Tot(i_h, j_h);
    
    
    % find optimal dvTot date indexes
    [minDV_opt, linearIdx] = min(deltavTot(:)); % Get the minimum value and its linear index
    [i_opt, j_opt] = ind2sub(size(deltavTot), linearIdx); % Convert linear index to (i, j) indices
    
    minDV1_opt = deltav1Tot(i_opt, j_opt);
    minDV2_opt = deltav2Tot(i_opt, j_opt);
    
    
    % % find Optimal dV2 (dv1max<maxDV1)
    % valid_mask_maxdv1 = (deltav1Tot < maxDV1);  % Logical(0 or 1) mask where deltav1Tot < maxDV1
    % valid_deltav2_maxdv1 = deltav2Tot(valid_mask_maxdv1);  % Extract valid deltav2Tot values in column vector 
    % 
    % if ~isempty(valid_deltav2_maxdv1)
    %     % Find the minimum deltav2Tot
    %     [minDV2_maxdv1, min_index] = min(valid_deltav2_maxdv1);
    % 
    %     % Find the corresponding (i, j) indices in the original array
    %     linear_indices_maxdv1 = find(valid_mask_maxdv1);  % Get linear indices of valid values
    %     min_valid_index_maxdv1 = linear_indices_maxdv1(min_index);  % Find the linear index of min deltav2Tot
    %     [i_maxdv1, j_maxdv1] = ind2sub(size(deltavTot), min_valid_index_maxdv1);  % Convert linear index to (i, j)
    % 
    %     % Retrieve the corresponding deltavTot value
    %     minDV_maxdv1 = deltavTot(i_maxdv1, j_maxdv1);
    % 
    %     % Retrieve the corresponding deltav1Tot value
    %     minDV1_maxdv1 = deltav1Tot(i_maxdv1, j_maxdv1);
    % else
    %     i_maxdv1 = 1;
    %     j_maxdv1 = 1;
    %     minDV_maxdv1 = 0;
    %     minDV1_maxdv1 = 0;
    %     minDV2_maxdv1 = 0;
    %     warning('_maxdv1 null')
    % end
    % 
    % desired_indexes = [i_h, j_h; ...
    %     i_opt, j_opt; ...
    %     i_maxdv1, j_maxdv1];
    % 
    % desired_minDV = [minDV_h, minDV1_h, minDV2_h;...
    %     minDV_opt, minDV1_opt, minDV2_opt;...
    %     minDV_maxdv1, minDV1_maxdv1, minDV2_maxdv1];

    desired_indexes = [i_h, j_h;  i_opt, j_opt];
    
    desired_minDV = [minDV_h, minDV1_h, minDV2_h;...
        minDV_opt, minDV1_opt, minDV2_opt];

end