function v1 = LMinETransfer(r1,r2,tm,mu)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LMinETransfer - Compute the initial velocity for a minimum-energy
%                 transfer
%
%   Inputs:
%       r1 - Initial position vector (km)
%       r2 - Final position vector (km)
%       tm - Time-of-flight mode indicator (typically 1 or -1)
%       mu - Gravitational parameter of the central body (km^3/s^2)
%
%   Outputs:
%       v1 - Initial velocity vector for the minimum-energy transfer (km/s)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % position norms
    R1 = norm(r1);
    R2 = norm(r2);

    % semi-minor axis c
    c = norm(r2-r1);

    % semilatus rectum ro
    costh = (dot(r1,r2))/(R1 * R2);
    pmin = ((R1*R2)/c) *(1-costh);
    
    F = 1 - (R2/pmin)*(1-costh);

    sinth = tm*sqrt(1-costh^2);
    G = (R2*R1)/sqrt(mu*pmin) * sinth;

    v1 = 1/G * (r2-F*r1);
end