% lambda_i and phi are in degrees
% R in meters

function [output] = ECEF_to_LLH(rx, ry ,rz)
radius_earth = 6378137; %from lecture slides
input = [rx, ry, rz];
R = norm(input);
lambda_i = asind(rz/R);
phi = atan2d(ry, rx); %atan2 takes into account quadrant
alt = R - radius_earth; 
output = [lambda_i, phi, alt];
end