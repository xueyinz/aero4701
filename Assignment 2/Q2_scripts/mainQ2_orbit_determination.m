%% 440305585
% AERO4701
%
% mainQ2_orbit_determination.m

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

% logical array for when satellite is visible
visible = pos_POLAR(ELEVATION,:) > min_elevation;

% times at which satellite is visible to ground station
t_visible = t_24hr;
t_visible(~visible) = NaN;

% portion of the satellite orbit visible to ground station
visible_POLAR_noisy = pos_POLAR_noisy;
visible_POLAR_noisy(:, ~visible) = NaN;
visible_ECI_noisy = polar2eci_vector(visible_POLAR_noisy, ground_LLH, t_visible);

%% Herrick-Gibbs technique for initial orbit determination

% estimate velocity at the second to the second last tracking measurements
initial_estimates = herrick_gibbs(visible_ECI_noisy, t_visible);

% logical array for when we have Herrick-Gibbs estimates
estimate_available = ~isnan(initial_estimates(1,:));

% time-stamps corresponding to each estimate
t_estimates = t_visible(estimate_available);

% use the first Herrick-Gibbs estimate for NLLS
first_estimate = initial_estimates(:, estimate_available);
first_estimate = first_estimate(:,1);
t_first_estimate = t_visible(estimate_available);
t_first_estimate = t_first_estimate(1);

% use the first Herrick-Gibbs estimate for position and velocity to get the orbital parameters
initial_parameters = posvel2orbital(first_estimate);
initial_satellite = get_satellite_struct(initial_parameters);
initial_satellite.OrbitalParameters = 'Initial (Herrick-Gibbs)';
initial_satellite.t0 = t_first_estimate;
initial_pos_ECI = orbit2ECI(initial_satellite, t_24hr - t_first_estimate);

%% Universal conic-section solution with NLLS method

% prepare input variables for nlls_orbit_determ script
observations = [t_visible(visible); visible_POLAR_noisy(:, visible)]';

% nlls method to refine the orbital parameters
refined_parameters = nlls_orbit_determ(observations, ground_LLH, first_estimate, t_first_estimate);
refined_satellite = get_satellite_struct(refined_parameters);
refined_satellite.OrbitalParameters = 'Refined (NLLS)';
refined_satellite.t0 = t_first_estimate;
refined_pos_ECI = orbit2ECI(refined_satellite, t_24hr - t_first_estimate);

%% Noise and errors

% noise profile
noise_profile = pos_POLAR_noisy - pos_POLAR;

% average satellite position errors
initial_error = rssq((initial_pos_ECI - pos_ECI), 1);
initial_error_avg = mean(initial_error);
refined_error = rssq((refined_pos_ECI - pos_ECI), 1);
refined_error_avg = mean(refined_error);