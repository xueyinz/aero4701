function LGCV = ecef2lgcv(ecef, groundstation_llh_geocentric)

x = ecef(1);
y = ecef(2);
z = ecef(3);
r = norm(ecef);
% phi = rad2deg(asin(-z/r));
% theta = rad2deg(atan2(y, x));
% POLAR = [r; phi; theta];

% latitude = groundstation(1);    % geocentric LLH
% longitude = groundstation(2);
% altitude = groundstation(3);

% convert groundstation geocentric LLH to ECEF
groundstation_ECEF = llh_geocentric2ecef(groundstation_llh_geocentric);
