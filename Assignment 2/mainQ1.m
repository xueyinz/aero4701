%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1.m

clearvars;
close all;

plotting = true;
save_figures = false;   % used for when saving graphs for the report

%% initialisation

fprintf('Performing calculations for Question 1...\n');
addpath('./Q1_data/', './Q1_useful_code/', './Q2_useful_code/', './coordinate_conversions/', './orbit_simulation/', './earth_plots/');
constants;              % load constants

%% Q1.A - satellite trajectory over a 12hr period

% load GPS ephemeris data
GPS_ephemeris = load_GPS_ephemeris('GPSsat_ephem.txt');
n_satellites = size(GPS_ephemeris, 1);
orbits.pos_ECI = zeros(3, size(t_12hr,2));
orbits = repmat(orbits, 1, n_satellites);

% for each satellite, get the trajectory over a 12hr period
for ii = 1:n_satellites
    
    % current satellite
    sat = GPS_ephemeris(ii);
    
    % mean anomaly: t = t1 - t0
    M = sat.M0*ones_t + sat.n.*(t_12hr - sat.t0);
    
    % eccentric anomaly
    E = mean2eccentric(M, sat.e);
    
    % true anomaly
    theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2));
    
    % radius
    r = sat.p*ones_t/(1 + sat.e.*cos(theta));
    
    % get ECI position vectors
    orbits(ii).pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);
    
end

%% Q1.B - satellite trajectory over the time period of the GPS pseudorange data

% load GPS pseudorange data for UAV
GPS_pseudorange = load_GPS_pseudorange('GPS_pseudorange.txt');   % [ time (s), sat #, pseudorange (s) ]
t_uav = [GPS_pseudorange.time];
ones_t_uav = ones(1, length(t_uav));
orbits_uav.pos_ECI = 0;
orbits_uav.pos_ECEF = 0;
orbits_uav = repmat(orbits_uav, 1, size(GPS_ephemeris, 1));

% for each satellite, get the trajectory at each unique time-stamp in the
% GPS pseudorange data
for ii = 1:n_satellites
    
    % current satellite
    sat = GPS_ephemeris(ii);
    
    % mean anomaly: t = t1 - t0
    M = sat.M0*ones_t_uav + sat.n.*(t_uav - sat.t0);
    
    % eccentric anomaly
    E = mean2eccentric(M, sat.e);
    
    % true anomaly
    theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2));
    
    % radius
    r = sat.p*ones_t_uav/(1 + sat.e.*cos(theta));
    
    % get ECI position vectors
    orbits_uav(ii).pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);
    
%     orbits_uav(ii).pos_ECEF = eci2ecef_vector(orbits_uav(ii).pos_ECI, sat.t0 - t_VE + t_uav);
        
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
    
%     if time_stamp == 210
%         pause;
%     end
    
    % current reading of pseudorange measurements
    reading = GPS_pseudorange(time_stamp);
    n_measurements = numel(reading.sat_num);
    
    % initial position estimate: centre of the Earth with no clock bias
    x_0 = [0; 0; 0; 0];
    
    % initial Jacobian matrix: one row for each satellite observation, 4
    % columns for each state parameter
    clear H;
    H = ones(n_measurements, 4);
    
    % residual errors
    clear delta_rho
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
            X_SV = orbits_uav(satellite).pos_ECI(:, time_stamp);
            
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
        delta_x = (transpose(H)*H)\transpose(H)*delta_rho;
        
        % update estimate for state x
        x_0 = x_0 + delta_x;
        
        % break out of NLLS update if error is sufficiently small
        if abs(delta_x) < error_threshold
            break;
        end
        
    end

    %get the degree of precisions
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
        fprintf('Bad rank\n');
        x_0 = [NaN; NaN; NaN; NaN];     % get rid of funky measurement
    end
    
    uav.pos_ECI(:,time_stamp) = x_0;
    
end

% convert uav coordinates into other frames
uav.pos_ECEF = eci2ecef_vector(uav.pos_ECI(XYZ,:), sat.t0 - t_VE + t_uav);
uav.pos_LLH = ecef2llh_geocentric_vector(uav.pos_ECEF(XYZ,:));
uav.pos_LGCV = ecef_ground2lgcv_vector(uav.pos_ECEF(XYZ,:), ground_LLH);
uav.pos_POLAR = cartesian2polar_vector(uav.pos_LGCV);
uav.pos_RANGE = uav.pos_POLAR(1,:);
uav.pos_AZ_deg = rad2deg(uav.pos_POLAR(2,:));           % convert from radians to degrees for azimuth
uav.pos_EL_deg = rad2deg(uav.pos_POLAR(3,:));           % convert from radians to degrees for elevation

% load actual UAV positions for comparison
counter = 0;
UAV_truth = load_UAV_position('UAVPosition.txt');
for tt = 1:numel(UAV_truth.time)
    if find(t_uav(t_uav == UAV_truth.time(tt)))
        counter = counter + 1;
        fprintf('%d\n', counter);
    end
end

%% plotting results for Question 1.A

if plotting == true
    
%     fprintf('Plotting 12hr orbit animation for Q1.A...\n');
% 
%     PlotEarthSphere;
%     view(3);
% 
%     % get all the plot handles for each satellite
%     for ii = 1:size(orbits,2)
%         plot3(orbits(ii).pos_ECI(1,:), orbits(ii).pos_ECI(2,:), orbits(ii).pos_ECI(3,:), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k');
%         orbits(ii).plot_handle = plot3(orbits(ii).pos_ECI(1,1), orbits(ii).pos_ECI(2,1), orbits(ii).pos_ECI(3,1), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
%         hold on;
%     end
% 
%     % loop through each time step and update position of each satellite
%     for ii = 1:size(t_12hr,2)
%         for jj = 1:size(orbits,2)
%             orbits(jj).plot_handle.XData = orbits(jj).pos_ECI(1,ii);
%             orbits(jj).plot_handle.YData = orbits(jj).pos_ECI(2,ii);
%             orbits(jj).plot_handle.ZData = orbits(jj).pos_ECI(3,ii);
%         end
%         drawnow;
%     end
% 
%     fprintf('Finished plotting animation for Q1.A.\nPlotting graphs for Q1.B...\n');

%% plotting results for Question 1.B

%     PlotEarthSphere;
%     view(3);
% 
%     for ii = 1:size(orbits_uav,2)
%         plot3(orbits_uav(ii).pos_ECI(1,:), orbits_uav(ii).pos_ECI(2,:), orbits_uav(ii).pos_ECI(3,:), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k');
%         orbits_uav(ii).plot_handle = plot3(orbits_uav(ii).pos_ECI(1,1), orbits_uav(ii).pos_ECI(2,1), orbits_uav(ii).pos_ECI(3,1), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
%         hold on;
%     end
% 
%     for ii = 1:size(t_uav,2)
%         for jj = 1:size(orbits_uav,2)
%             orbits_uav(jj).plot_handle.XData = orbits_uav(jj).pos_ECI(1,ii);
%             orbits_uav(jj).plot_handle.YData = orbits_uav(jj).pos_ECI(2,ii);
%             orbits_uav(jj).plot_handle.ZData = orbits_uav(jj).pos_ECI(3,ii);
%         end
%         drawnow;
%     end
% 
%     fprintf('Finished plotting graphs for Q1.B\n');

    figure;
    plot3(UAV_truth.LGCV(1,:), UAV_truth.LGCV(2,:), UAV_truth.LGCV(3,:),'.');
    hold on;
    plot3(uav.pos_LGCV(1,:), uav.pos_LGCV(2,:), uav.pos_LGCV(3,:), '.');
    
    title('UAV in LGCV coordinates')
    xlabel('X (m)')
    ylabel('Y (m)')
    zlabel('Z (m)')

end