%% 440305585
% AERO4701
%
% mainA1.m

clearvars;
close all;
% clc;

%% load data

addpath('./constants/', './tle/', './coordinate_conversions/', './orbit_simulation/', './Mitch_Bryson/');
globalConstants;
sat = extract_from_TLE('TLE_data.txt');

fprintf("Performing calculations for orbit simulation...\n");

%% simulate a 24h orbit starting from TLE Epoch

t = 0:dt:sec_per_day;                       % 24hr time vector; time since TLE Epoch in increments of dt [s]
ones_t = ones(1, length(t));                % array of ones of length t

M = sat.M0*ones_t + sat.n.*t;               % mean anomaly; note: t1 - t0 = t
E = mean2eccentric(M, sat.e);               % convert mean anomaly to eccentric anomaly by solving for E numerically
theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2)); % convert eccentric anomaly to true anomaly
r = sat.p*ones_t/(1 + sat.e.*cos(theta));   % calculate the radius
pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);     % get ECI position vectors

%% product a ground trace from the 24h orbit

pos_ECEF = eci2ecef_vector(pos_ECI, etime(sat.t0, t_VE) + t);   % convert ECI coordinates into ECEF coordinates
pos_LLH = ecef2llh_geocentric_vector(pos_ECEF);             % convert ECEF coordinates into LLH coordinates
pos_LL = rad2deg(pos_LLH(1:2, :));                         % convert latitude and longitude to degrees

%% plot the satellite orbit and ground trace

fprintf('Plotting satellite orbit  and ground trace...');

PlotEarthSphere;
plot3(pos_ECI(1,:),pos_ECI(2,:),pos_ECI(3,:), 'linewidth', 2, 'color', 'red');

% PlotEarthLatLong;
PlotEarthLatLong;
plot(pos_LL(2,:), pos_LL(1,:), 'r');

%% ground station perspective

% ground station
g