function [posErr, velErr] = odeError(dvData_selected, muSun) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% odeError - Compute propagation errors using an ODE solver
%
%   Inputs:
%       dvData_selected - Structure containing selected transfer entry
%                         with fields:
%                           date3_MJD  - Reentry departure date (MJD)
%                           date4_MJD  - Reentry arrival date (MJD)
%                           astID      - Asteroid identifier
%       muSun           - Gravitational parameter of the Sun (km^3/s^2)
%
%   Outputs:
%       posErr          - Norm of the position error (km) between the
%                         Lambert solution final position and the ODE
%                         solver final position
%       velErr          - Norm of the velocity error (km/s) between the
%                         Lambert solution
%                         final velocity and the ODE solver final velocity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    F = @(t,x) [x(4);               % dx/dt = Vx
                 x(5);              % dy/dt = Vy
                 x(6);              % dz/dt = Vz
                -muSun*x(1)/((x(1)^2+x(2)^2+x(3)^2)^(3/2));   % dVx/dt
                -muSun*x(2)/((x(1)^2+x(2)^2+x(3)^2)^(3/2));   % dVy/dt
                -muSun*x(3)/((x(1)^2+x(2)^2+x(3)^2)^(3/2) )]; % dVz/dt
    date3 = dvData_selected.date3_MJD;
    date4 = dvData_selected.date4_MJD;
    dT = (date4-date3)*86400;

    % seconds from departure to arrival
    tspan = [0, dT]; 
    
    % initial state (position and velocity from Lambert solver)
    astID = dvData_selected.astID;
    [kep, ~, ~] = ephNEO(date4, astID+11);
    p = kep(1) * (1 - kep(2)^2);
    carVec = kep2car(kep, muSun, p);
    r4 = carVec(1:3);
    [r3, ~] = EphSS_car(3, date3);
    [vSc_initial, vSc_final] = myLambert(r3, r4, dT, 1, muSun);
    
    y0 = [r3, vSc_initial];
    
    % ode
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [~,y] = ode45(F, tspan, y0, options);
    
    fpos = y(end, 1:3); % final postion
    vpos = y(end, 4:6); % final velocity
    
    % position
    posErr = norm(r4-fpos);
    velErr = norm(vSc_final-vpos);


end
