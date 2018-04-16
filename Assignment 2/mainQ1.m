%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1.m

clearvars;
close all;

save_figures = false;   % used for when saving graphs for the report

%% load data

addpath('./Q1_data/', './Q1_useful_code/', './Q2_useful_code/', './coordinate_conversions/', './orbit_simulation/', './earth_plots/');
global_constants;

% GPS ephemeris
GPS_ephemeris = load_GPS_ephemeris('GPSsat_ephem.txt');
orbit.pos_ECI = 0;
orbits = repmat(orbit, 1, size(GPS_ephemeris, 1));

% GPS pseudorange for UAV
GPS_pseudorange = load_GPS_pseudorange('GPS_pseudorange.txt');   % [ time (s), sat #, pseudorange (s) ]
t_uav = [GPS_pseudorange.time];
orbit_uav.pos_ECI = 0;
orbits_uav = repmat(orbit_uav, 1, size(GPS_pseudorange, 1));
% UAV.pos_ECI = 0;

%% Question 1.A

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
    orbits(ii).plot_handle = plot3(pos_ECI(1,1),pos_ECI(2,1),pos_ECI(3,1), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
end

%% Question 1.B

ones_t_uav = ones(1, length(t_uav));
for ii = 1:size(GPS_ephemeris, 1)
    sat = GPS_ephemeris(ii);
    M = sat.M0*ones_t_uav + sat.n.*(t_uav - sat.t0);    % mean anomaly; note: t1 - t0 = t
    E = mean2eccentric(M, sat.e);
%     theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2)) - pi/2; % convert eccentric anomaly to true anomaly
    theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2)); % convert eccentric anomaly to true anomaly
    r = sat.p*ones_t_uav/(1 + sat.e.*cos(theta));                   % calculate the radius
    pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);     % get ECI position vectors
    plot3(pos_ECI(1,:),pos_ECI(2,:),pos_ECI(3,:), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k');
    orbits_uav(ii).pos_ECI = pos_ECI;
    orbits_uav(ii).plot_handle = plot3(pos_ECI(1,1),pos_ECI(2,1),pos_ECI(3,1), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
end

%% plotting results

PlotEarthSphere;
view(3);

for jj = 1:size(t_12, 2)
    for kk = 1:size(orbits,2)
        orbits(kk).plot_handle.XData = orbits(kk).pos_ECI(1,jj);
        orbits(kk).plot_handle.YData = orbits(kk).pos_ECI(2,jj);
        orbits(kk).plot_handle.ZData = orbits(kk).pos_ECI(3,jj);
    end
    drawnow;
end