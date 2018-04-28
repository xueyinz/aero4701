%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1_B.m

%% Q1.B - multi-lateration of pseudorange data using non-linear least squares method

orbits_uav.pos_ECI = 0;
orbits_uav.pos_ECEF = 0;
orbits_uav = repmat(orbits_uav, 1, size(GPS_ephemeris, 1));
orbits_pos_LGCV = zeros(3, 2, n_satellites);
orbits_pos_POLAR = zeros(3, 2, n_satellites);

% for each satellite, get the trajectory at each unique time-stamp in the
% GPS pseudorange data
for ii = 1:n_satellites
    
    % current satellite
    sat = GPS_ephemeris(ii);
    
    % time since epoch
    t_since_epoch = t_uav - sat.t0;
    
    % get ECI position vectors
    orbits_uav(ii).pos_ECI = orbit2ECI(sat, t_since_epoch);
    
    % get ECEF position vectors
    orbits_uav(ii).pos_ECEF = eci2ecef_vector(orbits_uav(ii).pos_ECI, t_uav - t_VE);
    
    % get satellite position at start and end times (for plotting later)
    orbits_pos_LGCV(:,:,ii) = ecef_ground2lgcv_vector(orbits_uav(ii).pos_ECEF(:, [1,end]), ground_LLH);
    orbits_pos_POLAR(:,:,ii) = cartesian2polar_vector(orbits_pos_LGCV(:,:,ii));
        
end

%% estimation of UAV position using non-linear least squares of GPS pseudorange data

% initialisations
n_measurements = zeros(1, size(GPS_pseudorange,1));
GDOP = zeros(size(GPS_pseudorange,1),1);
PDOP = zeros(size(GPS_pseudorange,1),1);
HDOP = zeros(size(GPS_pseudorange,1),1);
VDOP = zeros(size(GPS_pseudorange,1),1);
TDOP = zeros(size(GPS_pseudorange,1),1);
uav.pos_ECEF = zeros(3, size(GPS_pseudorange,1));
uav.clock_bias = zeros(1, size(GPS_pseudorange,1));

% at each unique time-stamp in the pseudorange data, make a position
% estimate for the uav: x = [x, y, z, cb]
for time_stamp = 1:size(GPS_pseudorange,1)
    
    % current reading of pseudorange measurements
    reading = GPS_pseudorange(time_stamp);
    n_measurements(time_stamp) = numel(reading.sat_num);
    
    % initial position estimate: centre of the Earth with no clock bias
    x_0 = [0; 0; 0; 0];
    
    % initial Jacobian matrix: one row for each satellite observation, 4
    % columns for each state parameter
    H = ones(n_measurements(time_stamp), 4);
    
    % residual errors
    delta_rho = zeros(n_measurements(time_stamp), 1);
    
    % update estimation until convergence or we reach the maximum number of
    % iterations
    for jj = 1:max_iterations
        
        % update Jacobian matrix for each satellite measurement at the
        % current reading
        for measurement = 1:n_measurements(time_stamp)
            
            % number of current satellite of interest
            satellite = reading.sat_num(measurement);
            
            % ECEF position of relevant satellite at current time-stamp
            X_SV = orbits_uav(satellite).pos_ECEF(:, time_stamp);
            
            % temporary variable
            range = sqrt(sum((X_SV - x_0(XYZ)).^2));
            
            % estimate for range taking into account clock bias
            rho_estimate = range + x_0(CB);
            
            % update residual error
            delta_rho(measurement) = reading.range(measurement) - rho_estimate;
            
            % derivative of measurement rho w.r.t. state x
            % i.e. [drho/dx, drho/dy, drho/dz, drho/cb] where drho/cb = 1
            drho_dx = [-(X_SV - x_0(XYZ))/range; 1];
            
            % NLLS update of current measurement row in the Jacobian matrix
            H(measurement,:) = drho_dx';
            
        end
        
        % change in estimated state x
        delta_x = (transpose(H)*H) \ transpose(H) * delta_rho;
        
        % update estimate for state x
        x_0 = x_0 + delta_x;
        
        % break out of NLLS update if error is sufficiently small
        if abs(delta_x) < error_threshold
            break;
        end
        
    end

    % get degree of precision
    V = inv(transpose(H)*H);
    
    % dilutions of precision
    GDOP(time_stamp) = sqrt(V(1,1) + V(2,2) + V(3,3) + V(4,4));
    PDOP(time_stamp) = sqrt(V(1,1) + V(2,2) + V(3,3));
    HDOP(time_stamp) = sqrt(V(1,1) + V(2,2));
    VDOP(time_stamp) = sqrt(V(3,3));
    TDOP(time_stamp) = sqrt(V(4,4));
    
    % check rank of V
    if rank(V) ~= 4
        
        GDOP(time_stamp) = NaN;
        PDOP(time_stamp) = NaN;
        HDOP(time_stamp) = NaN;
        VDOP(time_stamp) = NaN;
        TDOP(time_stamp) = NaN;
        
        x_0 = [NaN; NaN; NaN; NaN];     % get rid of erroneous measurement
        
%         fprintf('Found erroneous measurement.\n');
        clc;
        
    end
    
    % store estimated state after removing bad value
    uav.pos_ECEF(:,time_stamp) = x_0(XYZ);
    uav.clock_bias(:, time_stamp) = x_0(CB);
    
end

%% coordinate transforms

% convert simulation truth of UAV coordinates into other frames
uav_truth.pos_POLAR = cartesian2polar_vector(uav_truth.pos_LGCV);

% convert estimated UAV coordinates into other frames
% uav.pos_LLH = ecef2llh_geocentric_vector(uav.pos_ECEF);
uav.pos_LGCV = ecef_ground2lgcv_vector(uav.pos_ECEF, ground_LLH);
uav.pos_POLAR = cartesian2polar_vector(uav.pos_LGCV);

%% removing outliers

% assume that satellite clock bias is negligible compared to user receiver
% clock bias. Remove estimations where clock bias is too large or too small
% because clock bias should be similar for the same GPS receiver
outliers = isoutlier(uav.clock_bias, 'grubbs');

uav.pos_ECEF_filtered = uav.pos_ECEF;
uav.pos_ECEF_filtered(:, outliers) = NaN;
uav.pos_ECEF_removed = uav.pos_ECEF;
uav.pos_ECEF_removed(:, ~outliers) = NaN;

% uav.pos_LLH_filtered = uav.pos_LLH;
% uav.pos_LLH_filtered(:, outliers) = NaN;
% uav.pos_LLH_removed = uav.pos_LLH;
% uav.pos_LLH_removed(:, ~outliers) = NaN;

uav.pos_LGCV_filtered = uav.pos_LGCV;
uav.pos_LGCV_filtered(:, outliers) = NaN;
uav.pos_LGCV_removed = uav.pos_LGCV;
uav.pos_LGCV_removed(:, ~outliers) = NaN;

uav.pos_POLAR_filtered = uav.pos_POLAR;
uav.pos_POLAR_filtered(:, outliers) = NaN;
uav.pos_POLAR_removed = uav.pos_POLAR;
uav.pos_POLAR_removed(:, ~outliers) = NaN;

n_measurements_filtered = n_measurements;
n_measurements_filtered(outliers) = NaN;
n_measurements_removed = n_measurements;
n_measurements_removed(~outliers) = NaN;

uav.clock_bias_filtered = uav.clock_bias;
uav.clock_bias_filtered(:, outliers) = NaN;
uav.clock_bias_removed = uav.clock_bias;
uav.clock_bias_removed(:, ~outliers) = NaN;

GDOP_filtered = GDOP;
PDOP_filtered = PDOP;
HDOP_filtered = HDOP;
VDOP_filtered = VDOP;
TDOP_filtered = TDOP;
GDOP_filtered(outliers) = NaN;
PDOP_filtered(outliers) = NaN;
HDOP_filtered(outliers) = NaN;
VDOP_filtered(outliers) = NaN;
TDOP_filtered(outliers) = NaN;
GDOP_removed = GDOP;
PDOP_removed = PDOP;
HDOP_removed = HDOP;
VDOP_removed = VDOP;
TDOP_removed = TDOP;
GDOP_removed(~outliers) = NaN;
PDOP_removed(~outliers) = NaN;
HDOP_removed(~outliers) = NaN;
VDOP_removed(~outliers) = NaN;
TDOP_removed(~outliers) = NaN;