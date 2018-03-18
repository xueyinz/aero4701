function LLH = ecef2llh_geocentric(pos_ECEF)

% convert Earth-Centred Earth-Fixed coordinate to geocentric LLH
r_Earth = 6378137;
x = pos_ECEF(1);
y = pos_ECEF(2);
z = pos_ECEF(3);
R = norm(pos_ECEF);
lambda = asin(z/R);     % latitude
phi = atan2(y, x);      % longitude
LLH = [rad2deg(lambda); rad2deg(phi); R - r_Earth];