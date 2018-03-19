%% 440305585
% AERO4701
%
% Convert Cartesian coordinates into polar coordinates
%
% inputs:   pos_CARTESIAN = [x; y; z] = Cartesian/LGCV coordinates [m, m, m]
% outputs:  pos_POLAR = [r, az, el] = polar coordinates [m, rad, rad]

function pos_POLAR = cartesian2polar(pos_CARTESIAN)

    if ((size(pos_CARTESIAN,1) ~= 3) || (size(pos_CARTESIAN,2) ~= 1))
        error('Check that dimensions of inputs to cartesian2polar match: ([x;y;z])')
    end

    x = pos_CARTESIAN(1);
    y = pos_CARTESIAN(2);
    z = pos_CARTESIAN(3);
    
    r = norm(pos_CARTESIAN);
    phi = atan2(y, x);          % azimuth
    theta = asin(-z/r);         % elevation
    
    pos_POLAR = [r; phi; theta];
    
end