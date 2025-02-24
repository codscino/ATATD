%% start
% Julian days conversion
t1 = date2mjd2000([2016, 3, 14, 12, 0, 0]); %earth time
t2 = date2mjd2000([2016, 10, 15, 12, 0, 0]); %mars time

% Positions
[r1,v1] = EphSS_car(3,t1);
[r2,v2] = EphSS_car(4,t2);

muSun = getAstroConstants('Sun','mu');

% velocity of the spacecraft
tm = 1;
vsc = LMinETransfer(r1,r2,tm,muSun);

%% compute final error
dT = 215*86400;
rSc_final = FGKepler_dt2(r1, vsc, dT, muSun) ;
dr = norm(rSc_final - r2);