%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1_B_plots.m

%% 3D plot

figure;
plot3(uav_truth.pos_LGCV(1,:), uav_truth.pos_LGCV(2,:), uav_truth.pos_LGCV(3,:), 'k-');
hold on;
plot3(uav.pos_LGCV_filtered(1,:), uav.pos_LGCV_filtered(2,:), uav.pos_LGCV_filtered(3,:), '.');
plot3(uav.pos_LGCV_removed(1,:), uav.pos_LGCV_removed(2,:), uav.pos_LGCV_removed(3,:), 'r.');
grid on;
% axis equal;
title('3D plot of UAV position w.r.t. to ground station (LGCV)');
xlabel('X (north) (m)');
ylabel('Y (east) (m)');
zlabel('Z (down) (m)');
legend('Simulation truth', 'Estimations from GPS pseudorange', 'Filtered out estimations', 'Location', 'northeast');
if save_figures == true
    saveas(gcf, 'Q1B_3D_plot.png');
    view(2);
    saveas(gcf, 'Q1B_3D_plot_xy.png');
    view(3);
end

%% 2D overhead plot

figure;
[plotx_truth, ploty_truth] = overheadplotcoords(uav_truth.pos_POLAR(AZIMUTH,:), uav_truth.pos_POLAR(ELEVATION,:));
[plotx, ploty] = overheadplotcoords(uav.pos_POLAR_filtered(AZIMUTH,:), uav.pos_POLAR_filtered(ELEVATION,:));
[plotx_removed, ploty_removed] = overheadplotcoords(uav.pos_POLAR_removed(AZIMUTH,:), uav.pos_POLAR_removed(ELEVATION,:));
plot(plotx_truth, ploty_truth, 'k-');
hold on;
plot(plotx, ploty, '.');
plot(plotx_removed, ploty_removed, 'r.');
grid on;
axis equal;
title('2D overhead plot of UAV w.r.t. ground station (polar)');
xlabel('Azimuth (degrees)');
ylabel('Elevation (degrees)');
legend('Simulation truth', 'Estimations from GPS pseudorange', 'Filtered out estimations', 'Location', 'eastoutside');
if save_figures == true
    saveas(gcf, 'Q1B_2D_overhead_plot.png');
end

%% altitude vs. time

figure;
plot(uav_truth.time, uav_truth.pos_LGCV(ALTITUDE,:), 'k');
hold on;
plot(t_uav, uav.pos_LGCV_filtered(ALTITUDE,:), '.');
plot(t_uav, uav.pos_LGCV_removed(ALTITUDE,:), 'r.');
grid on;
title('Altitude w.r.t. ground station vs. time (lgcv)');
xlabel('Epoch time (s)');
ylabel('Altitude (m)');
legend('Simulation truth', 'Estimations from GPS pseudorange', 'Filtered out estimations', 'Location', 'southwest');
xlim([min(t_uav) max(t_uav)]);
if save_figures == true
    saveas(gcf, 'Q1B_altitude_vs_time.png');
end

%% DOP

figure;
plot(t_uav, GDOP_filtered, 'b');
hold on
plot(t_uav, PDOP_filtered, 'g');
plot(t_uav, HDOP_filtered, 'r');
plot(t_uav, VDOP_filtered, 'm');
plot(t_uav, TDOP_filtered, 'k');
plot(t_uav, GDOP_removed, 'b.');
plot(t_uav, PDOP_removed, 'g.');
plot(t_uav, HDOP_removed, 'r.');
plot(t_uav, VDOP_removed, 'm.');
plot(t_uav, TDOP_removed, 'k.');
grid on;
xlim([min(t_uav) max(t_uav)]);
title('Dilution of precision vs. time');
xlabel('Epoch time (s)');
ylabel('Dilution of precision');
legend('GDOP', 'PDOP', 'HDOP', 'VDOP', 'TDOP', 'GDOP filtered out', 'PDOP filtered out', 'HDOP filtered out', 'VDOP filtered out', 'TDOP filtered out', 'Location', 'northwest');
if save_figures == true
    saveas(gcf, 'Q1B_DOP.png');
end

%% number of satellite used

figure;
bar(t_uav, n_measurements_filtered);
hold on;
% bar(t_uav, n_measurements_removed, 'FaceColor', [grey grey grey], 'EdgeColor', [grey grey grey]);
bar(t_uav, n_measurements_removed, 'FaceColor', 'r', 'EdgeColor', 'r');
title('Number of satellites used at each time-stamp');
xlabel('Epoch time (s) (discrete)');
ylabel('Number of satellites used');
legend('Readings used in estimation', 'Readings removed from estimation', 'Location', 'southeast');
if save_figures == true
    saveas(gcf, 'Q1B_n_satellites_used.png');
end

%% user clock bias

figure;
plot(t_uav, uav.clock_bias_filtered, '.');
hold on;
plot(t_uav, uav.clock_bias_removed, 'r.');
grid on;
xlim([min(t_uav) max(t_uav)]);
title('User clock bias');
xlabel('Epoch time (s)');
ylabel('User clock bias (m)');
legend('Estimations from GPS pseudorange', 'Filtered out estimations', 'Location', 'southwest');
if save_figures == true
    saveas(gcf, 'Q1B_user_clock_bias.png');
end

%% azimuth-elevation of satellites visible to GPS receiver w.r.t. ground station (polar)

figure;
visible_sat_start = GPS_pseudorange(1).sat_num;
visible_sat_end = GPS_pseudorange(end).sat_num;
theta_start = orbits_pos_POLAR(AZIMUTH, 1, visible_sat_start);
theta_end = orbits_pos_POLAR(AZIMUTH, end, visible_sat_end);
r_start = orbits_pos_POLAR(ELEVATION, 1, visible_sat_start);
r_end = orbits_pos_POLAR(ELEVATION, end, visible_sat_end);
polarscatter(theta_start, 90 - rad2deg(r_start));
rticks([0 30 60 90]);
rlim([0 90]);
hold on;
polarscatter(theta_end, 90 - rad2deg(r_end));
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.ThetaTickLabels = {'Az = N','30deg','60','E','120','150','S','210','240','W','300','330'};
ax.ThetaTickLabels = {'Az = N','30','60','E','120','150','S','210','240','W','300','330'};
ax.RTickLabels = {'El = 90deg','60','30','0'};
title('Azimuth-elevation of visible satellites (polar w.r.t. ground station)');
legend('Satellites visible at the start','Satellites visible at the end');
if save_figures == true
    saveas(gcf, 'Q1B_azimuth-elevation.png');
end