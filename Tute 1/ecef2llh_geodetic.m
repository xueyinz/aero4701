%% 440305585
% AERO4701
%
% Convert Earth-Centred Earth-Fixed coordinates to geodetic latitude,
% longitude, height coordinates
%
% inputs:   pos_ECEF = [x; y; z] = ECEF coordinates [m, m, m]
% outputs:  pos_LLH = [lat; long; alt] = LLH geodetic coordinates [rad, rad, m]

function LLH = ecef2llh_geodetic(pos_ECEF)

    if ((size(pos_ECEF,1) ~= 3) || (size(pos_ECEF,2) ~= 1))
        error('Check that dimensions of inputs to ecef2llh_geodetic match: ([x;y;z])')
    end
    
    global a;
    global e;

    x = pos_ECEF(1);
    y = pos_ECEF(2);
    z = pos_ECEF(3);

    % initial values
    h = 0;
    N = a;
    p = sqrt(x^2 + y^2);
    lambda = 100;

    % iterate to solve for latitude, longitude, altitude
    for i = 1:100
        sinlambda = z/(N*(1 - e^2) + h);
        lambda = atan((z + (e^2)*N*sinlambda)/p);
        N = a/sqrt(1 - (e^2)*(sin(lambda)^2));
        h = p/cos(lambda) - N;
    end

    LLH = [lambda; atan2(y,x); h];

end