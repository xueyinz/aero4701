%% 440305585
% AERO4701
%
% Convert Earth-Centred Earth-Fixed coordinates to geocentric latitude,
% longitude, height coordinates
%
% inputs:   pos_ECEF = [x; y; z] = ECEF coordinates [m, m, m]
% outputs:  pos_LLH = [lat; long; alt] = LLH geocentric coordinates [rad, rad, m]

function pos_LLH = ecef2llh_geocentric(pos_ECEF)

    if ((size(pos_ECEF,1) ~= 3) || (size(pos_ECEF,2) ~= 1))
        error('Check that dimensions of inputs to eci2ecef match: ([x;y;z])')
    end

    global r_Earth;
    
    x = pos_ECEF(1);
    y = pos_ECEF(2);
    z = pos_ECEF(3);
    
    R = norm(pos_ECEF);   % altitude
    lambda = asin(z/R);             % latitude
    phi = atan2(y, x);              % longitude
    h = R - r_Earth;
    
    pos_LLH = [lambda; phi; h];

end