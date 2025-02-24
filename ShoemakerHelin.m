function [dv, astClass] = ShoemakerHelin(a, e, i_deg)
    % asteroidDeltaV calculates the delta-v for rendezvous with an asteroid
    % and classifies the asteroid.
    %
    %   [dv, astClass] = ShoemakerHelin(a, e, i_deg)
    %
    % Inputs:
    %   a     - semimajor axis [AU]
    %   e     - eccentricity (dimensionless)
    %   i_deg - inclination [deg]
    %
    % Outputs:
    %   dv      - estimated delta-v [km/s]
    %   astClass- asteroid classification:
    %             1 if Amors (1.017 AU < q < 1.3 AU)
    %             2 if Apollos (q < 1.017 AU)
    %             3 if Atens (a < 1.0 AU and Q > 0.983 AU)
    
    % Convert inclination from degrees to radians
    i = i_deg * pi/180;   % [rad]
    
    % Compute perihelion (q) and aphelion (Q) distances [AU]
    q = a * (1 - e);      % [AU]
    Q = a * (1 + e);      % [AU]
    
    % Classify the asteroid and compute the corresponding formulas:
    if a < 1.0
        % Atens: a < 1.0 AU and Q > 0.983 AU
        if Q > 0.983
            astClass = 3; % 3 for Atens
            % Use Aten formulas
            ut2 = 2 - 2*cos(i/2)*sqrt(2*Q - Q^2);               
            uc2 = 3/Q - 1 - (2/Q)*sqrt(2 - Q);                    
            ur2 = 3/Q - 1/a - (2/Q)*cos(i/2)*sqrt(a*(1 - e^2)/Q);   
        else
            error('For Aten asteroids, aphelion Q must be > 0.983 AU.');
        end
    elseif a > 1.0
        % For a > 1.0 AU, distinguish Apollos vs. Amors based on q.
        if q < 1.017
            % Apollos: q < 1.017 AU
            astClass = 2; % 2 for Apollos
            ut2 = 3 - 2/(Q + 1) - 2*cos(i/2)*sqrt(2*Q/(Q+1));      
            uc2 = 3/Q - 2/(Q+1) - (2/Q)*sqrt(2/(Q+1));              
            ur2 = 3/Q - 1/a - (2/Q)*cos(i/2)*sqrt((a/Q)*(1 - e^2));   
        elseif (q > 1.017) && (q < 1.3)
            % Amors: 1.017 AU < q < 1.3 AU
            astClass = 1; % 1 for Amors
            ut2 = 3 - 2/(Q + 1) - 2*cos(i/2)*sqrt(2*Q/(Q+1));      
            uc2 = 3/Q - 2/(Q+1) - (2/Q)*cos(i/2)*sqrt(2/(Q+1));     
            ur2 = 3/Q - 1/a - (2/Q)*sqrt(a*(1 - e^2)/Q);             
        else
            error('For a > 1.0 AU, q must be either < 1.017 AU (Apollo) or between 1.017 and 1.3 AU (Amor).');
        end
    else
        error('Asteroid classification not defined for the given parameters.');
    end
    
    % Constants and normalized velocities
    v_earth = 29.784;       % Earth's orbital velocity [km/s]
    U0 = 3.074 / v_earth;     % Normalized geostationary velocity (dimensionless)
    S  = sqrt(2) * U0;        % Normalized escape velocity from geostationary orbit (dimensionless)
    
    % Impulse for leaving geostationary orbit.
    ul = sqrt(ut2 + S^2) - U0;  
    
    % Impulse for rendezvousing at the asteroid.
    ur = sqrt(uc2 - 2*sqrt(ur2*uc2)*cos(i/2) + ur2);  
    
    % Figure of merit from Shoemaker and Helin.
    F = ul + ur;  % dimensionless
    
    % Delta-v calculation (scaling converts F to km/s)
    dv = 30 * F + 0.5;  % [km/s]
end