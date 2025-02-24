function [vSc_initial, vSc_final] = myLambert(r1,r2,dTarget,tm,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% myLambert - Solve Lambert's problem using iterative refinement
%
%   Inputs:
%       r1      - Initial position vector (km)
%       r2      - Final position vector (km)
%       dTarget - Target time of flight (seconds)
%       tm      - Time-of-flight mode indicator (typically 1 or -1)
%       mu      - Gravitational parameter of the central body (km^3/s^2)
%
%   Outputs:
%       vSc_initial - Computed initial spacecraft velocity vector (km/s)
%       vSc_final   - Computed final spacecraft velocity vector (km/s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    tol = 10^(-3); %1 metre tolerance
    N = 20; %total number of iterations for the for loop
    maxNumIter = 25; % max number of iterations for the while loop

    [~, ~, dtmin] = MinETransfer(r1,r2,tm,mu); %dtmin from the min energy transfer(231 days)
    dT = dTarget; %target time of 215 days

    vSc_initial = LMinETransfer(r1,r2,tm,mu);
    
    for i = 1:N %i index of the for loop
        lambda = i/N;
        dt = lambda*dT + (1-lambda)*dtmin;
    
        Error = 10;
        numIter = 0; % iterations while loops
        
        while (Error>tol) && (numIter<maxNumIter)
            numIter = numIter + 1;

            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% checks and warning %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % a check
            a = mu / ((2*mu)/norm(r1) - norm(vSc_initial)^2);
            if a < 0
                %warning("DOH! a is negative (a = %f), the method diverges", a);
                vSc_initial = nan(1,3);
                vSc_final = nan(1,3);
                return
            end

            % iterations check
            if numIter > maxNumIter
                % warning("DOH! max number of iterations reached, the method diverges (numIter = %d)", numIter);
                vSc_initial = nan(1,3);
                vSc_final = nan(1,3);
                return
            end


            %%%%%%%%%%%%%%%%%%%%%%%
            %%% finish warnings %%%
            %%%%%%%%%%%%%%%%%%%%%%%

            [rSc_final, ~] = FGKepler_dt2(r1, vSc_initial, dt, mu);
            Smat = STM_Lambert(r1, vSc_initial, dt, mu);
            dr_t2 = r2 - rSc_final;
            dv_t1 = inv(Smat) * dr_t2';
            vSc_initial = vSc_initial + dv_t1';
            Error = norm(dr_t2);
        end
    end

    
    
    % recheck if a is negative before calculating vSc_final 

    a_final = mu / ((2*mu)/norm(r1) - norm(vSc_initial)^2);

    if a_final < 0
         % warning("DOH! a is negative (a = %f), the method diverges", a);
         vSc_initial = nan(1,3);
         vSc_final = nan(1,3);
         return
    end

    [~, vSc_final] = FGKepler_dt2(r1, vSc_initial, dt, mu);
    

end