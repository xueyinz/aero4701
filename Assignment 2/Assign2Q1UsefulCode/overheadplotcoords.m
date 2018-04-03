function [plotx, ploty] = overheadplotcoords(azi, elev)

% OVERHEADPLOTCOORDS - computes the MATLAB x and y coordinates for and
% overhead azimuth and elevation plot (input in radians)

plotx = -(0.5*pi - elev)*(180/pi)*sin(azi);
ploty = (0.5*pi - elev)*(180/pi)*cos(azi);

% end of overheadplotcoords

