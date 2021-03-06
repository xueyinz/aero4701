%% 440305585
% AERO4701
%
% mainQ2_plots.m

%% orbital parameters

print_orbital_parameters(sat, initial_satellite, refined_satellite);


%%  24hr period orbit trajectories

PlotEarthSphere;
view(3);

plot3(pos_ECI(1,:), pos_ECI(2,:), pos_ECI(3,:), 'r');
grid on;
hold on;
plot3(initial_pos_ECI(1,:), initial_pos_ECI(2,:), initial_pos_ECI(3,:), 'g');
plot3(refined_pos_ECI(1,:), refined_pos_ECI(2,:), refined_pos_ECI(3,:), 'b');

xlim([min([pos_ECI(1,:) initial_pos_ECI(1,:) refined_pos_ECI(1,:)]) max([pos_ECI(1,:) initial_pos_ECI(1,:) refined_pos_ECI(1,:)])])
ylim([min([pos_ECI(2,:) initial_pos_ECI(2,:) refined_pos_ECI(2,:)]) max([pos_ECI(2,:) initial_pos_ECI(2,:) refined_pos_ECI(2,:)])])
zlim([min([pos_ECI(3,:) initial_pos_ECI(3,:) refined_pos_ECI(3,:)]) max([pos_ECI(3,:) initial_pos_ECI(3,:) refined_pos_ECI(3,:)])])
legend('Earth', 'Reference trajectory', 'Initial orbit determination', 'Orbit refinement', 'Location', 'eastoutside');

if save_figures == true
    saveas(gcf, 'Q2_24hr_trajectories.png');
end

%% noise profile

figure;
subplot(3,1,1);
plot(t_24hr, noise_profile(1,:));
grid on;
title('Noise profile');
ylabel('Range error (m)');

subplot(3,1,2);
plot(t_24hr, rad2deg(noise_profile(2,:)));
grid on;
ylabel('Azimuth error (deg)');

subplot(3,1,3);
plot(t_24hr, rad2deg(noise_profile(3,:)));
grid on;
ylabel('Elevation error (deg)');
xlabel('Time (s)');

if save_figures == true
    saveas(gcf, 'Q2_noise_profile.png');
end

%% tracking measurements

figure;

subplot(3,1,1);
plot(t_24hr, visible_POLAR_noisy(1,:), 'r.');
hold on;
plot(t_24hr, pos_POLAR(1,:), '-');
grid on;
title('Ground station tracking measurements');
ylabel('Range (m)');

subplot(3,1,2);
plot(t_24hr, rad2deg(visible_POLAR_noisy(2,:)), 'r.');
hold on;
plot(t_24hr, rad2deg(pos_POLAR(2,:)), '-');
grid on;
ylabel('Azimuth (deg)');

subplot(3,1,3);
plot(t_24hr, rad2deg(visible_POLAR_noisy(3,:)), 'r.');
hold on;
plot(t_24hr, rad2deg(pos_POLAR_noisy(3,:)), '-');
grid on;
ylabel('Elevation (deg)');
xlabel('Time (s)');
legend('Noisy observations', 'Actual orbit', 'Location', 'northeast');

if save_figures == true
    saveas(gcf, 'Q2_tracking_measurements.png');
end

%% satellite position error

fprintf('Average satellite position error for initial orbit determination: %.2fm\n', initial_error_avg);
fprintf('Average satellite position error for refined orbit determination: %.2fm\n', refined_error_avg);

figure;
plot(t_24hr, initial_error);
grid on;
hold on;
plot(t_first_estimate, initial_error(find(t_24hr == t_first_estimate)), 'or');
title('Satellite position error: initial orbit determination');
xlabel('Time (s)');
ylabel('Position error (m)');
legend('24 hr orbit','First estimate only', 'Location', 'southeast');

if save_figures == true
    saveas(gcf, 'Q2_error_initial.png');
end

figure;
plot(t_24hr, refined_error);
grid on;
title('Satellite position error: refined orbit determination');
xlabel('Time (s)');
ylabel('Position  error (m)');

if save_figures == true
    saveas(gcf, 'Q2_error_refined.png');
end

figure;
plot(t_24hr, initial_error);
grid on;
hold on;
plot(t_24hr, refined_error);
plot(t_first_estimate, initial_error(find(t_24hr == t_first_estimate)), 'or');
title('Satellite position error');
xlabel('Time (s)');
ylabel('Position error (m)');
legend('Initial estimate', 'Refined estimate', 'First estimate only', 'Location', 'southeast');

if save_figures == true
    saveas(gcf, 'Q2_error_both.png');
end