function dvData2 = Leg2ConditionsPar(dvData1, date4_max, step3, ...
    ToF2, RE, muEarth, muSun)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Leg2ConditionsPar - Compute reentry leg maneuver parameters for transfers
%
%   Inputs:
%       dvData1   - Structure array containing first-stage transfer entries
%       date4_max - Maximum allowable reentry arrival date (MJD)
%       step3     - Step size (days) for constructing the reentry
%                   departure grid
%       ToF2      - Array of time-of-flight values for the reentry
%                   leg (days)
%       RE        - Earth's radius (km)
%       muEarth   - Gravitational parameter of Earth (km^3/s^2)
%       muSun     - Gravitational parameter of the Sun (km^3/s^2)
%
%   Outputs:
%       dvData2   - Structure array containing valid reentry transfer
%                   entries with fields:
%                     date1_MJD  - First-stage departure date (MJD)
%                     date2_MJD  - First-stage arrival date (MJD)
%                     date3_MJD  - Reentry departure date (MJD)
%                     date4_MJD  - Reentry arrival date (MJD)
%                     dv1        - Delta-v for first-stage departure (km/s)
%                     dv2        - Delta-v for first-stage arrival (km/s)
%                     dv3        - Delta-v for reentry departure (km/s)
%                     dv4        - Delta-v for reentry arrival (km/s)
%                     dvLeg1     - Total delta-v for first stage (km/s)
%                     dvLeg2     - Total delta-v for reentry leg (km/s)
%                     dv234      - Sum of dv2, dv3, and dv4 (km/s)
%                     dvTot      - Total delta-v for the complete
%                                  transfer (km/s)
%                     ToF        - Total time-of-flight from first-stage
%                                  departure to reentry arrival (days)
%                     layover    - Layover time between first-stage arrival
%                                  and reentry departure (days)
%                     astID      - Asteroid identifier
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get the asteroid ID (assuming dvData1 is non-empty)
    astID = dvData1(1).astID;
    eph_astID = astID + 11;
    
    nEntries = length(dvData1);
    % Preallocate a cell array to hold sub-results from each dvData1 entry
    dvData2Cells = cell(nEntries, 1);
    
    parfor k = 1:nEntries
        subResults = [];  % Initialize an empty structure array for this iteration
        % Retrieve first-stage parameters for the k-th entry
        dv1_val = dvData1(k).deltav1;
        dv2_val = dvData1(k).deltav2;
        dvLeg1_val = dvData1(k).deltavLeg1;
        date1   = dvData1(k).date1_MJD;
        date2   = dvData1(k).date2_MJD;
        
        % Build the reentry departure grid for this first-stage transfer:
        date3_start = date2 + 60;   % minimum 2 months after arrival
        date3_end   = date2 + 180;    % maximum 6 months after arrival
        grid3 = date3_start:step3:date3_end;
        L3 = length(grid3);
        L4 = length(ToF2);
        
        % Pre-calculate the orbital parameter "p" at the start of grid3
        [kep, ~, ~] = ephNEO(grid3(1), eph_astID);
        p = kep(1) * (1 - kep(2)^2);
        
        % Loop over the reentry departure grid and second-stage
        % time-of-flight grid.
        for i3 = 1:L3
            for j4 = 1:L4
                date3 = grid3(i3);
                date4 = date3 + ToF2(j4);
                if date4 > date4_max
                    continue;
                end
                
                % Get asteroid state at reentry departure.
                [kep_reentry, ~, ~] = ephNEO(date3, eph_astID);
                carVec = kep2car(kep_reentry, muSun, p);
                r3 = carVec(1:3);
                v3 = carVec(4:6);
                
                % Get Earth state at reentry arrival.
                [r4, v4] = EphSS_car(3, date4);
                dT = ToF2(j4) * 86400; % Transfer time in seconds
                
                % --- Lambert solution with tm = 1 ---
                % [vSc_initial, vSc_final, ~, ~] = izzoLambert(r3, r4, ...
                %     dT/86400, 0, muSun);
                [vSc_initial, vSc_final] = myLambert(r3, r4, dT, 1, muSun);
                dv3_short = norm(vSc_initial - v3);
                v_inf_Earth = norm(v4 - vSc_final);
                Rf = RE + 400;
                Vpe = sqrt(2*(muEarth/Rf + v_inf_Earth^2/2));
                Vce = sqrt(muEarth/Rf);
                dv4_short = Vpe - Vce;
                deltavLeg2_short = dv3_short + dv4_short;
                
                % --- Lambert solution with tm = -1 ---
                % [vSc_initialM, vSc_finalM, ~, ~] = izzoLambert(r3, r4, ...
                %     -dT/86400, 0, muSun);
                [vSc_initialM, vSc_finalM] = myLambert(r3, r4,dT,-1,muSun);
                dv3_long = norm(vSc_initialM - v3);
                v_inf_EarthM = norm(v4 - vSc_finalM);
                VpeM = sqrt(2*(muEarth/Rf + v_inf_EarthM^2/2));
                dv4_long = VpeM - Vce;
                deltavLeg2_long = dv3_long + dv4_long;
                
                % Choose the solution with the minimum delta-v for reentry.
                if deltavLeg2_short < deltavLeg2_long
                    dv3 = dv3_short;
                    dv4 = dv4_short;
                    dvLeg2 = deltavLeg2_short;
                else
                    dv3 = dv3_long;
                    dv4 = dv4_long;
                    dvLeg2 = deltavLeg2_long;
                end
                
                dv234 = dv2_val + dv3 + dv4;
                
                % Only keep entries with dv3 < 1.
                if dv3 < 0.8
                    newEntry = struct(...
                        'date1_MJD', date1, ...
                        'date2_MJD', date2, ...
                        'date3_MJD', date3, ...
                        'date4_MJD', date4, ...
                        'dv1', dv1_val, ...
                        'dv2', dv2_val, ...
                        'dv3', dv3, ...
                        'dv4', dv4, ...
                        'dvLeg1', dvLeg1_val,...
                        'dvLeg2', dvLeg2, ...
                        'dv234', dv234, ...
                        'dvTot',  dv1_val+dv234, ...
                        'ToF', date4-date1, ...
                        'layover', date3-date2, ...
                        'astID', astID);
                    subResults = [subResults; newEntry]; 
                end
            end
        end
        
        dvData2Cells{k} = subResults;
    end
    
    % Concatenate the results from all iterations into one structure array.
    dvData2 = vertcat(dvData2Cells{:});
end