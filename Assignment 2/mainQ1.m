%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1.m

clearvars;
close all;

save_figures = false;   % used for when saving graphs for the report

%% load data

addpath('./Assign2Q1Data/', './Assign2Q1UsefulCode/', './Assign2Q2UsefulCode/', './coordinate_conversions/', './orbit_simulation/', './Mitch_Bryson/');
global_constants;

GPS_ephemeris = load_GPS_ephemeris('GPSsat_ephem.txt');
orbit.pos_ECI = 0;
orbits = repmat(orbit, 1, size(GPS_ephemeris, 1));
% orbits = ones(size(GPS_ephemeris,1), length(t_12));

PlotEarthSphere;
view(3);

for ii = 1:size(GPS_ephemeris, 1)
    sat = GPS_ephemeris(ii);
    M = sat.M0*ones_t + sat.n.*(t_12 - sat.t0);    % mean anomaly; note: t1 - t0 = t
    E = mean2eccentric(M, sat.e);
%     theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2)) - pi/2; % convert eccentric anomaly to true anomaly
    theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2)); % convert eccentric anomaly to true anomaly
    r = sat.p*ones_t/(1 + sat.e.*cos(theta));                   % calculate the radius
    pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);     % get ECI position vectors
    plot3(pos_ECI(1,:),pos_ECI(2,:),pos_ECI(3,:), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k');
    orbits(ii).pos_ECI = pos_ECI;
    orbits(ii).plot_handle = plot3(pos_ECI(1,1),pos_ECI(2,1),pos_ECI(3,1), 'o', 'MarkerEdgeColor','r', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
end

for jj = 1:size(t_12, 2)
    for kk = 1:size(orbits,2)
        orbits(kk).plot_handle.XData = orbits(kk).pos_ECI(1,jj);
        orbits(kk).plot_handle.YData = orbits(kk).pos_ECI(2,jj);
        orbits(kk).plot_handle.ZData = orbits(kk).pos_ECI(3,jj);
    end
    drawnow;
end