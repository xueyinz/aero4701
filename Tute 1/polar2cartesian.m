function CARTESIAN = polar2cartesian(polar)

r = polar(1);               % range
phi = deg2rad(polar(2));    % elevation
theta = deg2rad(polar(3));  % azimuth
x = r*cos(theta)*cos(phi);
y = r*cos(theta)*sin(phi);
z = -r*sin(theta);
CARTESIAN = [x; y; z];