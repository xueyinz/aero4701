%% 440305585
% AERO4701
% Assignment 1
% 
% Run me!
%
% mainA1.m

clearvars;
close all;
% clc;

save_figures = false;   % used for when saving graphs for the report

%% load data

addpath('./constants/', './tle/', './coordinate_conversions/', './orbit_simulation/', './Mitch_Bryson/');
globalConstants;
sat = extract_from_TLE('TLE_data.txt');                     % satellite parameters
start_epoch = sprintf("%02d/%02d/%4d %02d:%02d:%02d UT", sat.day, sat.month, sat.year, sat.hour, sat.min, sat.sec);
end_epoch = sprintf("%02d/%02d/%4d %02d:%02d:%02d UT", sat.day + 1, sat.month, sat.year, sat.hour, sat.min, sat.sec);
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

% plot 3D orbit simulation in ECI
PlotEarthSphere;
max_axis = max(max(abs(pos_ECI)));
axis([-max_axis max_axis -max_axis max_axis -max_axis max_axis]);
view(3);
plot3(pos_ECI(1,:),pos_ECI(2,:),pos_ECI(3,:), 'LineWidth', 1, 'Color', 'r');
legend('Earth', 'Satellite orbit', 'Location', 'eastoutside');
if save_figures
    saveas(gcf,'orbit_3D_ECI.png');
end
p1 = plot3(pos_ECI(1,1),pos_ECI(2,1),pos_ECI(3,1), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 5);
legend('Earth', 'Satellite orbit', 'Satellite instantaneous position');

% plot 3D orbit simulation in ECEF
PlotEarthSphere_ECEF;
max_axis = max(max(abs(pos_ECEF)));
axis([-max_axis max_axis -max_axis max_axis -max_axis max_axis]);
view(3);
plot3(pos_ECEF(1,:),pos_ECEF(2,:),pos_ECEF(3,:), 'LineWidth', 1, 'Color', 'r');
legend('Earth', 'Satellite orbit', 'Location', 'eastoutside');
if save_figures
    saveas(gcf,'orbit_3D_ECEF.png');
end

% plot ground trace
PlotEarthLatLong;
temp_ii = 1;
for ii = 1:(length(t)-1)    % plot without the horizontal lines when longitude crosses over from +ve to -ve or vice versa
    if (abs(pos_LL_deg(2,ii)) >  90) && (xor((pos_LL_deg(2,ii) < 0), (pos_LL_deg(2,ii+1) < 0)))
        p2 = plot(pos_LL_deg(2,temp_ii:ii), pos_LL_deg(1,temp_ii:ii), 'linewidth', 1, 'color', 'r');
        temp_ii = ii+1;
    end
end
axis image;
p3 = plot(ground_LLH_deg(2), ground_LLH_deg(1), 'd', 'MarkerFaceColor','g', 'MarkerEdgeColor','g');
legend([p2, p3], {'Satellite ground trace', 'Ground station'}, 'Location', 'eastoutside');
if save_figures
    saveas(gcf,'orbit_ground_trace.png');
end
p4 = plot(pos_LL_deg(2,1), pos_LL_deg(1,1), 'o', 'MarkerFaceColor','y', 'MarkerEdgeColor','y', 'MarkerSize', 5);
legend([p2, p3, p4], {'Satellite ground trace', 'Ground station', 'Satellite instantaneous position'}, 'Location', 'eastoutside');

%% ground station perspective

fprintf("Performing calculations for ground station perspective...\n");

ground_LLH = [deg2rad(ground_LLH_deg(1:2)); 0];             % ground station
pos_LGCV = ecef_ground2lgcv_vector(pos_ECEF, ground_LLH);   % convert ECEF coordinates into LGCV coordinates relative to ground station
pos_POLAR = cartesian2polar_vector(pos_LGCV);               % convert LGCV coordinates into polar coordinates relative to ground station
pos_RANGE = pos_POLAR(1,:);
pos_AZ_deg = rad2deg(pos_POLAR(2,:));                       % convert from radians to degrees for azimuth
pos_EL_deg = rad2deg(pos_POLAR(3,:));                       % convert from radians to degrees for elevation

max_el = max(pos_EL_deg);
t_sat_in_view = length(pos_EL_deg(pos_EL_deg >= 0)) * dt;   % sum of time during which satellite is in view, i.e. elevation >= 0 (horizon)
t_ratio_in_view = t_sat_in_view/sec_per_day;                % ratio of the 24h orbit that the satellite is visible to ground station

fprintf("\t* Ground station: %.2f latitude [deg], %.2f longitude [deg], %.2f altitude [m]\n", ground_LLH_deg(1), ground_LLH_deg(2), ground_LLH_deg(3));
fprintf("\t* Maximum elevation of satellite: %.2f [deg]\n", max_el);
fprintf("\t* Percentage time satellite is visible: %.2f%%\n\n", t_ratio_in_view*100);

%% plot the ground station perspective of the satellite

fprintf("Plotting satellite orbit wrt the ground station...\n");

t_visible = t(pos_EL_deg >= 0);                             % satellite visible when angle of elevation is positive
t_invisible = t(pos_EL_deg < 0);

% range
figure;
subplot(3,1,1);
axis([t(1) t(end) min(pos_RANGE) max(pos_RANGE)]);
scatter(t_visible, pos_RANGE(pos_EL_deg >= 0), scatter_size, 'g');
grid on;
hold on;
scatter(t_invisible, pos_RANGE(pos_EL_deg < 0), scatter_size, 'r');
title('Satellite range vs time since Epoch');
xlabel('Time since Epoch (seconds)');
ylabel('Satellite range (m)');
legend('Satellite visible', 'Satellite invisible', 'Location', 'eastoutside');

% azimuth
subplot(3,1,2);
axis([t(1) t(end) -180 180]);
scatter(t_visible, pos_AZ_deg(pos_EL_deg >= 0), scatter_size, 'g');
grid on;
hold on;
scatter(t_invisible, pos_AZ_deg(pos_EL_deg < 0), scatter_size, 'r');
title('Satellite azimuth vs time since Epoch');
xlabel('Time since Epoch (seconds)');
ylabel('Satellite azimuth (deg)');
legend('Satellite visible', 'Satellite invisible', 'Location', 'eastoutside');

% elevation
subplot(3,1,3);
axis([t(1) t(end) -90 90]);
scatter(t_visible, pos_EL_deg(pos_EL_deg >= 0), scatter_size, 'g');
grid on;
hold on;
scatter(t_invisible, pos_EL_deg(pos_EL_deg < 0), scatter_size, 'r');
title('Satellite elevation vs time since Epoch');
xlabel('Time since Epoch (seconds)');
ylabel('Satellite elevation (deg)');
legend('Satellite visible', 'Satellite invisible', 'Location', 'eastoutside');
if save_figures
    saveas(gcf,'ground_station.png');
end

% azimuth vs elevation
figure;
scatter(pos_AZ_deg(pos_EL_deg >= 0), pos_EL_deg(pos_EL_deg >= 0), scatter_size, 'g');
axis([-180 180 -90 90]);
grid on;
hold on;
scatter(pos_AZ_deg(pos_EL_deg < 0), pos_EL_deg(pos_EL_deg < 0), scatter_size, 'r');
title('Satellite azimuth vs satellite elevation');
xlabel('Satellite azimuth (deg)');
ylabel('Satellite elevation (deg)');
legend('Satellite visible', 'Satellite invisible', 'Location', 'eastoutside');
if save_figures
    saveas(gcf,'ground_az_vs_el.png');
end

%% run animations

% plot the moving satellite animation
for jj = 2:length(pos_ECI(1,:))
    p1.XData = pos_ECI(1,jj);
    p1.YData = pos_ECI(2,jj);
    p1.ZData = pos_ECI(3,jj);
    p4.XData = pos_LL_deg(2,jj);
    p4.YData = pos_LL_deg(1,jj);
    drawnow
end