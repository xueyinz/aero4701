%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1_B.m

orbits_uav.pos_ECI = 0;
orbits_uav.pos_ECEF = 0;
orbits_uav = repmat(orbits_uav, 1, size(GPS_ephemeris, 1));

% for each satellite, get the trajectory at each unique time-stamp in the
% GPS pseudorange data
for ii = 1:n_satellites
    
    % current satellite
    sat = GPS_ephemeris(ii);
    
    % mean anomaly: t = t1 - t0
    M = sat.M0 + sat.n.*(t_uav - sat.t0);
    
    % eccentric anomaly
    E = mean2eccentric(M, sat.e);
    
    % true anomaly
    theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2));
    
    % radius
    r = sat.p./(1 + sat.e.*cos(theta));
    
    % get ECI position vectors
    orbits_uav(ii).pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);
    
    % get ECEF position vectors
    orbits_uav(ii).pos_ECEF = eci2ecef_vector(orbits_uav(ii).pos_ECI, t_uav - t_VE);
        
end

%% Q1.B - multi-lateration of pseudorange data using non-linear least squares method

% initialisations
GDOP = zeros(size(GPS_pseudorange,1),1);
PDOP = zeros(size(GPS_pseudorange,1),1);
HDOP = zeros(size(GPS_pseudorange,1),1);
VDOP = zeros(size(GPS_pseudorange,1),1);
TDOP = zeros(size(GPS_pseudorange,1),1);
uav.pos_ECEF = zeros(4, size(GPS_pseudorange,1));

% at each unique time-stamp in the pseudorange data, make a position
% estimate for the uav: x = [x, y, z, cb]
for time_stamp = 1:size(GPS_pseudorange,1)
    
    % current reading of pseudorange measurements
    reading = GPS_pseudorange(time_stamp);
    n_measurements = numel(reading.sat_num);
    
    % initial position estimate: centre of the Earth with no clock bias
    x_0 = [0; 0; 0; 0];
    
    % initial Jacobian matrix: one row for each satellite observation, 4
    % columns for each state parameter
    H = ones(n_measurements, 4);
    
    % residual errors
    delta_rho = zeros(n_measurements, 1);
    
    % update estimation until convergence or we reach the maximum number of
    % iterations
    for jj = 1:max_iterations
        
        % update Jacobian matrix for each satellite measurement at the
        % current reading
        for measurement = 1:n_measurements
            
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
        
        clc;
        x_0 = [NaN; NaN; NaN; NaN];     % get rid of erroneous measurement
        
    end
    
    uav.pos_ECEF(:,time_stamp) = x_0;
    
end

% convert uav coordinates into other frames
uav.pos_LLH = ecef2llh_geocentric_vector(uav.pos_ECEF(XYZ,:));
uav.pos_LGCV = ecef_ground2lgcv_vector(uav.pos_ECEF(XYZ,:), ground_LLH);
uav.pos_POLAR = cartesian2polar_vector(uav.pos_LGCV);
uav.pos_RANGE = uav.pos_POLAR(1,:);
uav.pos_AZ_deg = rad2deg(uav.pos_POLAR(2,:));           % convert from radians to degrees for azimuth
uav.pos_EL_deg = rad2deg(uav.pos_POLAR(3,:));           % convert from radians to degrees for elevation
