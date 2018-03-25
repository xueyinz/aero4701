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
start_epoch = sprintf("%02d/%02d/%4d %02d:%02d:%02d", sat.day, sat.month, sat.year, sat.hour, sat.min, sat.sec);
end_epoch = sprintf("%02d/%02d/%4d %02d:%02d:%02d", sat.day + 1, sat.month, sat.year, sat.hour, sat.min, sat.sec);
fprintf("\nOrbit of %s over a 24hr period from %s to %s\n\n", sat.name, start_epoch, end_epoch);

%% simulate a 24h orbit starting from TLE Epoch

fprintf("Performing calculations for orbit simulation...\n");

t = 0:dt:sec_per_day;                                       % 24hr time vector; time since TLE Epoch in increments of dt [s]
ones_t = ones(1, length(t));                                % array of ones of length t

M = sat.M0*ones_t + sat.n.*t;                               % mean anomaly; note: t1 - t0 = t
E = mean2eccentric(M, sat.e);                               % convert mean anomaly to eccentric anomaly by solving for E numerically
theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2)); % convert eccentric anomaly to true anomaly
r = sat.p*ones_t/(1 + sat.e.*cos(theta));                   % calculate the radius

pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);     % get ECI position vectors
pos_ECEF = eci2ecef_vector(pos_ECI, etime(sat.t0, t_VE) + t);   % convert ECI coordinates into ECEF coordinates
pos_LLH = ecef2llh_geocentric_vector(pos_ECEF);             % convert ECEF coordinates into LLH coordinates
pos_LL_deg = rad2deg(pos_LLH(1:2, :));                      % convert latitude and longitude to degrees

period_min = floor(sat.period/60);
period_sec = floor(sat.period - period_min*60);
fprintf("\t* Orbital period: %dmin %dsec\n\n", period_min, period_sec);

%% plot the satellite orbit and ground trace

fprintf('Plotting satellite orbit and ground trace...\n\n');

% plot 3D orbit simulation
PlotEarthSphere;
plot3(pos_ECI(1,:),pos_ECI(2,:),pos_ECI(3,:), 'linewidth', 2, 'color', 'red');

% plot ground trace
PlotEarthLatLong;
plot(pos_LL_deg(2,:), pos_LL_deg(1,:), 'r');

%% ground station perspective

fprintf("Performing calculations for ground station perspective...\n");

ground_LLH_deg = [-23.7; 133.87; 0];                        % ground station in LLH [deg, deg, m]
ground_LLH = [deg2rad(ground_LLH_deg(1:2)); 0];             % ground station
pos_LGCV = ecef_ground2lgcv_vector(pos_ECEF, ground_LLH);   % convert ECEF coordinates into LGCV coordinates relative to ground station
pos_POLAR = cartesian2polar_vector(pos_LGCV);               % convert LGCV coordinates into polar coordinates relative to ground station
pos_AZ_deg = rad2deg(pos_POLAR(2,:));                       % convert from radians to degrees for azimuth
pos_EL_deg = rad2deg(pos_POLAR(3,:));                       % convert from radians to degrees for elevation

max_el = max(pos_EL_deg);
t_sat_in_view = length(pos_EL_deg(pos_EL_deg >= 0)) * dt;   % sum of time during which satellite is in view, i.e. elevation >= 0 (horizon)
t_ratio_in_view = t_sat_in_view/sec_per_day;                % ratio of the 24h orbit that the satellite is visible to ground station

fprintf("\t* Ground station: %.2f latitude [deg], %.2f longitude [deg], %.2f altitude [m]\n", ground_LLH_deg(1), ground_LLH_deg(2), ground_LLH_deg(3));
fprintf("\t* Maximum elevation of satellite: %.2f [deg]\n", max_el);
fprintf("\t* Percentage time satellite is visible: %.2f%%\n\n", t_ratio_in_view*100);

%% plot the ground station perspective of the satellite

fprintf("Plotting satellite orbit wrt the ground station\n");

