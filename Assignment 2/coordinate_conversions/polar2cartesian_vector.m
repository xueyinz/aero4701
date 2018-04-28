%% 440305585
% AERO4701
%
% Convert polar coordinates into Cartesian coordinates
%
% inputs:   pos_POLAR = [r, az, el] = polar coordinates [m, rad, rad]
% outputs:	pos_CARTESIAN = [x; y; z] = Cartesian/LGCV coordinates [m, m, m]

function pos_CARTESIAN = polar2cartesian_vector(pos_POLAR)

    r = pos_POLAR(1,:);
    phi = pos_POLAR(2,:);             % elevation
    theta = pos_POLAR(3,:);           % azimuth
    
    x = r.*cos(theta).*cos(phi);
    y = r.*cos(theta).*sin(phi);
    z = -r.*sin(theta);
    
    pos_CARTESIAN = [x; y; z];

end