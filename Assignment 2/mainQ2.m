%% 440305585
% AERO4701
% Assignment 2
%
% mainQ2.m

clearvars;
close all;

plotting = false;
save_figures = false;   % used for when saving graphs for the report

%% initialisation

addpath('./Q2_scripts/', './coordinate_conversions/', './orbit_simulation/');

% load constants
Q2_constants;

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

% get ECEF position vectors - assume that t_VE = epoch (i.e. time of last
% vernal equinox alignment = t0)
pos_ECEF = eci2ecef_vector(pos_ECI, t_24hr);

% get LGCV position vectors
pos_LGCV = ecef_ground2lgcv_vector(pos_ECEF, ground_LLH);

% get polar position vectors
pos_POLAR = cartesian2polar_vector(pos_LGCV);

% portion of the satellite orbit visible to ground station
visible_POLAR = pos_POLAR(:, (pos_POLAR(ELEVATION,:) > min_elevation));

% times at which satellite is visible to ground station
t_visible = t_24hr(pos_POLAR(ELEVATION,:) > min_elevation);


error = -deviation + rand()*2*deviation;


% plot results
if plotting == true
    
    % plot the 24hr orbit
    PlotEarthSphere;
    view(3);
    plot3(pos_ECI(1,:), pos_ECI(2,:), pos_ECI(3,:));
    
end
