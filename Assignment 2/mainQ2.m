%% 440305585
% AERO4701
% Assignment 2
%
% mainQ2.m

clearvars;
close all;

plotting = true;
save_figures = true;   % used for when saving graphs for the report

%% initialisation

addpath('./Q2_scripts/', './coordinate_conversions/', './orbit_simulation/');

% load constants
constants;

% mean anomaly
M = sat.M0 + sat.n.*(t_24hr);

% eccentric anomaly
E = mean2eccentric(M, sat.e);

% true anomaly
theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2));

% radius
r = sat.p./(1 + sat.e.*cos(theta));

% get ECI position vectors
pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);

% plot the 24hr orbit
PlotEarthSphere;
view(3);
plot3(pos_ECI(1,:), pos_ECI(2,:), pos_ECI(3,:));