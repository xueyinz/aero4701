%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1_B_plots.m

%% 3D plot

figure;
plot3(uav_truth.pos_LGCV(1,:), uav_truth.pos_LGCV(2,:), uav_truth.pos_LGCV(3,:), '.');
hold on;
plot3(uav.pos_LGCV(1,:), uav.pos_LGCV(2,:), uav.pos_LGCV(3,:), 'r.');
grid on;
% axis equal;
title('3D plot of UAV position w.r.t. to ground station (LGCV)');
xlabel('X (north) (m)');
ylabel('Y (east) (m)');
zlabel('Z (down) (m)');
legend('Simulation truth', 'Estimation from GPS pseudorange', 'Location', 'northeast');
if save_figures == true
    saveas(gcf, '3D_plot.png');
    view(2);
    saveas(gcf, '3D_plot_xy.png');
    view(3);
end

%% 2D overhead plot

figure;
[plotx_truth, ploty_truth] = overheadplotcoords(uav_truth.pos_POLAR(AZIMUTH,:), uav_truth.pos_POLAR(ELEVATION,:));
[plotx, ploty] = overheadplotcoords(uav.pos_POLAR(AZIMUTH,:), uav.pos_POLAR(ELEVATION,:));
plot(plotx_truth, ploty_truth, '.');
hold on;
plot(plotx, ploty, 'ro');
grid on;
axis equal;
title('2D overhead plot of UAV w.r.t. ground station (polar)');
xlabel('Azimuth (degrees)');
ylabel('Elevation (degrees)');
legend('Simulation truth', 'Estimation from GPS pseudorange', 'Location', 'eastoutside');
if save_figures == true
    saveas(gcf, '2D_overhead_plot.png');
end

%% altitude vs. time

figure;
plot(t_uav, uav.pos_LLH(ALTITUDE,:));
title('Altitude w.r.t. ground station vs. time (LLH)');
xlabel('Epoch time (s)');
ylabel('Altitude (m)');
grid on;
xlim([min(t_uav) max(t_uav)]);
if save_figures == true
    saveas(gcf, 'altitude_vs_time.png');
end

%% DOP

figure;
plot(t_uav, GDOP);
hold on
plot(t_uav, PDOP);
plot(t_uav, HDOP);
plot(t_uav, VDOP);
plot(t_uav, TDOP);
grid on;
xlim([min(t_uav) max(t_uav)]);
title('Dilution of precision vs. time');
xlabel('Epoch time (s)');
ylabel('Dilution of precision');
legend('GDOP', 'PDOP', 'HDOP', 'VDOP', 'TDOP', 'Location', 'northwest');
if save_figures == true
    saveas(gcf, 'DOP.png');
end

%% number of satellite used

figure;
bar(t_uav, n_measurements_filtered);
hold on;
grey = 0.6;
bar(t_uav, n_measurements_removed, 'FaceColor', [grey grey grey], 'EdgeColor', [grey grey grey]);
title('Number of satellites used at each time-stamp');
xlabel('Epoch time (s) (discrete)');
ylabel('Number of satellites used');
legend('Readings used in estimation', 'Readings removed from estimation', 'Location', 'southeast');
if save_figures == true
    saveas(gcf, 'n_satellites_used.png');
end