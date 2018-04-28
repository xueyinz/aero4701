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

%% Get clean satellite orbit

% get ECI position vectors
pos_ECI = orbit2ECI(sat, t_24hr);

% get ECEF position vectors - assume that t_VE = epoch (i.e. time of last
% vernal equinox alignment = t0)
pos_ECEF = eci2ecef_vector(pos_ECI, t_24hr);

% get LGCV position vectors
pos_LGCV = ecef_ground2lgcv_vector(pos_ECEF, ground_LLH);

% get polar position vectors
pos_POLAR = cartesian2polar_vector(pos_LGCV);

%% Add random error to the observations

% noisy observations in polar coordinates
pos_POLAR_noisy = pos_POLAR;
pos_POLAR_noisy(AZ_EL,:) = add_random_error(pos_POLAR_noisy(AZ_EL,:), angle_error);
pos_POLAR_noisy(RANGE,:) = add_random_error(pos_POLAR_noisy(RANGE,:), range_error);

% noisy observations in ECI frame
pos_ECI_noisy = polar2eci_vector(pos_POLAR_noisy, ground_LLH, t_24hr);

%% Ground station perspective of satellite

% times at which satellite is visible to ground station
t_visible = t_24hr;
t_visible(pos_POLAR(ELEVATION,:) < min_elevation) = NaN;

% portion of the satellite orbit visible to ground station
visible_POLAR_noisy = pos_POLAR_noisy;
visible_POLAR_noisy(:, (pos_POLAR(ELEVATION,:) < min_elevation)) = NaN;
visible_ECI_noisy = polar2eci_vector(visible_POLAR_noisy, ground_LLH, t_visible);

%% Herrick-Gibbs technique for initial orbit determination

% estimate velocity at the second to the second last tracking measurements
v_estimates = herrick_gibbs(pos_ECI_noisy, t_visible);

%% Universal conic-section solution with NLLS method

% prepare input variables for nlls_orbit_determ script
observations = [t_visible(~isnan(t_visible)); visible_POLAR_noisy(:, ~isnan(t_visible))];
ground_ECEF = llh_geocentric2ecef(ground_LLH);


% plot results
if plotting == true
    
    % plot the 24hr orbit
    PlotEarthSphere;
    view(3);
    plot3(pos_ECI(1,:), pos_ECI(2,:), pos_ECI(3,:));
    
end
