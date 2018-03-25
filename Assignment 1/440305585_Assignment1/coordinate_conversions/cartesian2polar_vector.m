%% 440305585
% AERO4701
%
% Convert Cartesian coordinates into polar coordinates
%
% inputs:   pos_CARTESIAN = [x; y; z] = Cartesian/LGCV coordinates [m, m, m]
% outputs:  pos_POLAR = [r, az, el] = polar coordinates [m, rad, rad]

function pos_POLAR = cartesian2polar_vector(pos_CARTESIAN)

    x = pos_CARTESIAN(1,:);
    y = pos_CARTESIAN(2,:);
    z = pos_CARTESIAN(3,:);
    
%     r = norm(pos_CARTESIAN);  % returns scalar for vectors and matrices
    r = sqrt(sum(pos_CARTESIAN.^2, 1));
    phi = atan2(y, x);          % azimuth
    theta = asin(-z./r);        % elevation
    
    pos_POLAR = [r; phi; theta];
    
end