function dvData = Leg1ConditionsPar(astID, grid1, ToF, v_earth, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leg1ConditionsPar - Compute transfer maneuver parameters for Leg 1
%
%   Inputs:
%       astID   - Asteroid identifier (used to index ephemeris data)
%       grid1   - Array of departure dates (MJD)
%       ToF     - Array of time-of-flight values (days)
%       v_earth - Earth's velocity at departure (km/s)
%       mu      - Gravitational parameter of the central body (km^3/s^2)
%
%   Outputs:
%       dvData  - Structure array containing valid transfer entries with
%                 fields:
%                   date1_MJD  - Departure date (MJD)
%                   date2_MJD  - Arrival date (MJD)
%                   deltav1    - Departure leg delta-v (km/s)
%                   deltav2    - Arrival leg delta-v (km/s)
%                   deltavLeg1 - Total delta-v for the transfer (km/s)
%                   astID      - Asteroid identifier
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L1 = length(grid1);
    L2 = length(ToF);
    
    % Pre-calculate the asteroid's orbital parameter "p" 
    eph_astID = astID + 11;
    [kep, ~, ~] = ephNEO(grid1(1), eph_astID );
    p = kep(1) * (1 - kep(2)^2);
    
    % Total number of iterations across the grid.
    N = L1 * L2;
    results = cell(N, 1);
    
    % Use parfor to loop over all (i1, j2) pairs.
    parfor idx = 1:N
        % Map idx to corresponding i1 and j2
        i1 = floor((idx-1) / L2) + 1;
        j2 = mod(idx-1, L2) + 1;
        
        % Earth state at departure
        date1 = grid1(i1);
        [r1, v1] = EphSS_car(3, date1);
        
        % Arrival time computed as departure + time-of-flight
        date2 = date1 + ToF(j2);
        
        % Asteroid state at arrival
        [kep_arr, ~, ~] = ephNEO(date2, eph_astID);
        carVec = kep2car(kep_arr, mu, p);
        r2 = carVec(1:3);
        v2 = carVec(4:6);
        
        % Compute transfer time in seconds
        dT = ToF(j2) * 86400;
        
        % --- Lambert solution with tm = 1 ---
        % [vSc_initial, vSc_final, ~, ~] = izzoLambert(r1, r2, ...
        %     dT/86400, 0, mu);
        [vSc_initial, vSc_final] = myLambert(r1, r2, dT, 1, mu);
        deltav1_short = norm(vSc_initial - v1);
        deltav2_short = norm(v2 - vSc_final);
        
        % --- Lambert solution with tm = -1 ---
        % [vSc_initialM, vSc_finalM, ~, ~] = izzoLambert(r1, r2, ...
        %     -dT/86400, 0, mu);
        [vSc_initialM, vSc_finalM] = myLambert(r1, r2, dT, -1, mu);
        deltav1_long = norm(vSc_initialM - v1);
        deltav2_long = norm(v2 - vSc_finalM);
        
        % Total delta-v for both solutions
        deltavLeg1_short = deltav1_short + deltav2_short;
        deltavLeg1_long  = deltav1_long  + deltav2_long;
        
        % Choose the minimum of the two solutions.
        if deltavLeg1_short < deltavLeg1_long
            dv1 = deltav1_short;
            dv2 = deltav2_short;
            dvLeg1 = deltavLeg1_short;
        else
            dv1 = deltav1_long;
            dv2 = deltav2_long;
            dvLeg1 = deltavLeg1_long;
        end
        
        % Check validity conditions:
        %   - dvLeg1 must be a number (not NaN)
        %   - dv1 must be less than 1.5 and dv2 less than 0.5
        %   - Additionally, check the chosen Lambert solution's
        %     vSc_final vs. v_earth.
        valid = false;
        if ~isnan(dvLeg1) && (dv1 < 1.5) && (dv2 < 0.5)
            if deltav1_short <= deltav1_long
                valid = abs(norm(vSc_final) - v_earth) < 1.5;
            else
                valid = abs(norm(vSc_finalM) - v_earth) < 1.5;
            end
        end
        
        % If the transfer is valid, store the result.
        if valid
            newEntry = struct('date1_MJD', date1, ...
                              'date2_MJD', date2, ...
                              'deltav1',   dv1, ...
                              'deltav2',   dv2, ...
                              'deltavLeg1', dvLeg1, ...
                              'astID',     astID);
            results{idx} = newEntry;
        else
            results{idx} = [];
        end
    end
    
    % Concatenate non-empty cells to form the output structure array.
    dvData = [results{:}];
end