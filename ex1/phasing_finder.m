function [hohmann_date, Thohmann, phasing_target] = phasing_finder(start_date, end_date)

    %%% Input %%%
    % start_date - Start date in MJD2000
    % end_date - End date in MJD2000
    
    %%% Output %%%
    % hohmann_date - MJD2000 date of optimal Hohmann transfer
    % Thohmann - Semi-period of the Hohmann transfer orbit (days)
    
    % Constants
    muSun = getAstroConstants('Sun', 'mu');  % Gravitational parameter of the Sun
    
    % aEarth = getAstroConstants('Earth', 'Sma'); % Semi-major axis of Earth
    aMars = getAstroConstants('Mars', 'Sma');   % Semi-major axis of Mars
    % 
    %% Compute Synodic Period and Hohmann Transfer Period
    %nEarth = sqrt(muSun / aEarth^3); % Mean motion of Earth [rad/s]
    nMars = sqrt(muSun / aMars^3);   % Mean motion of Mars [rad/s] 
    % Ts = (2 * pi) / abs(nMars - nEarth); % Synodic period [s]
    
    % Hohmann transfer semi-period
    rEarth = 150e6; % Earth-Sun mean distance [km]
    rMars = 1.524 * rEarth; % Mars-Sun mean distance [km]
    aH = (rMars + rEarth) / 2; % Hohmann transfer semi-major axis
    TH = 0.5 * 2 * pi * sqrt(aH^3 / muSun); % Hohmann transfer semi-period
    
    Thohmann = TH / 86400; % Convert seconds to days

    %% Compute Required Phase Angle
    % phasing angle at departure 
    phasing_target = mod(180 - (nMars * TH * 180/pi), 360); %Expected 44 degrees
    
    %% Iterate to Find Transfer Window
    step = 1; % Step in days
    tolerance = 0.5; % Tolerance in degrees
    t = start_date;

    while t <= end_date
        % Get planetary positions
        [rEarthVec, ~] = EphSS_car(3, t); % Earth position vector
        [rMarsVec, ~] = EphSS_car(4, t);  % Mars position vector
        
        % Compute heliocentric longitudes
        lambdaEarth = atan2d(rEarthVec(2), rEarthVec(1)); % Earth longitude
        lambdaMars = atan2d(rMarsVec(2), rMarsVec(1));   % Mars longitude
        
        % Compute actual phase angle (modulo 360)
        theta_actual = mod(lambdaMars - lambdaEarth, 360);
        
        % Check if within tolerance
        if abs(theta_actual - phasing_target) < tolerance
            hohmann_date = t;
            return;
        end
        
        % Increment day
        t = t + step;
    end

    % If no valid date is found
    error('No valid Hohmann transfer date found within the given range.');
end