% STM(state transition matrix) lambert
function Smat = STM_Lambert(r0, v0, dt, mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STM_Lambert - Compute the state transition matrix (STM) for Lambert's problem
%
%   Inputs:
%       r0  - Initial position vector (km)
%       v0  - Initial velocity vector (km/s)
%       dt  - Time interval (seconds)
%       mu  - Gravitational parameter of the central body (km^3/s^2)
%
%   Outputs:
%       Smat - State transition matrix (3x3) that relates changes in the 
%              initial state to changes in the final state.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R0 = norm(r0);
    V0 = norm(v0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% this is basically FGKepler_dt function %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    a = mu / ((2*mu)/R0 - V0^2);

    % claudio check
    % if isnan(a)
    %     Smat = nan(3,3);
    %     warning('smat is nan')
    %     return
    % end

    n = sqrt(mu/a^3);
    dM = n*dt;
    sigma0 = dot(r0,v0) / sqrt(mu);

    % find dE with fzero
    fun = @(dE) -dM + dE - (1- (R0/a))*sin(dE) - (sigma0/sqrt(a))*(cos(dE)-1);
    x = fzero(fun,dM);
    dE = x;

    F = 1 - (a/R0)*(1-cos(dE));
    G = dt + sqrt(a^3/mu)*(sin(dE)-dE);

    rf = F*r0 + G*v0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% end of FGKepler_dt function %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Rf = norm(rf);
    dF = - (sqrt(mu*a)/(Rf*R0)) * sin(dE);
    dG = 1 - (a/Rf) *(1-cos(dE));
    vf = dF*r0 + dG*v0;

    C = a*sqrt(a^3/mu) * (3*sin(dE) - (2+cos(dE)) *dE) -a*dt*(1-cos(dE));
    deltar = rf-r0;
    deltav = vf-v0;
    
    % sensitivity matrix
    Smat = (R0/mu) * (1-F)*(deltar'*v0 - deltav'*r0) + ...
    (C/mu)*vf'*v0 + G*eye(3);
end