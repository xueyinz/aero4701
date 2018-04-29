%% 440305585
% AERO4701
%
% mainQ2_secondary.m

%% experiment with sampling time

% save results from timestep = 100
t_24hr_100 = t_24hr;
initial_error_100 = initial_error;
refined_error_100 = refined_error;

% new results
t_24hr = 0:50:sec_per_day;
mainQ2_orbit_determination;

% plot together
figure;
plot(t_24hr, initial_error);
grid on;
hold on;
plot(t_24hr_100, initial_error_100);
title('Satellite position error: initial orbit determination');
xlabel('Time (s)');
ylabel('Position error (m)');
legend('50s time-step','100s time-step', 'Location', 'southeast');

if save_figures == true
    saveas(gcf, 'Q2_error_initial_timestep.png');
end

figure;
plot(t_24hr, refined_error);
grid on;
hold on;
plot(t_24hr_100, refined_error_100);
legend('50s time-step','100s time-step', 'Location', 'southeast');
title('Satellite position error: refined orbit determination');
xlabel('Time (s)');
ylabel('Position  error (m)');

if save_figures == true
    saveas(gcf, 'Q2_error_refined_timestep.png');
end

%% experiment with noise

range_error = range_error*10;
angle_error = angle_error*10;
mainQ2_orbit_determination;

% plot together
figure;
plot(t_24hr_100, initial_error_100);
grid on;
hold on;
plot(t_24hr, initial_error);
title('Satellite position error: initial orbit determination');
xlabel('Time (s)');
ylabel('Position error (m)');
legend('Standard noise','10x noise', 'Location', 'southeast');

if save_figures == true
    saveas(gcf, 'Q2_error_initial_noise.png');
end

figure;
plot(t_24hr_100, refined_error_100);
grid on;
hold on;
plot(t_24hr, refined_error);
legend('Standard noise','10x noise', 'Location', 'southeast');
title('Satellite position error: refined orbit determination');
xlabel('Time (s)');
ylabel('Position  error (m)');

if save_figures == true
    saveas(gcf, 'Q2_error_refined_noise.png');
end

% %% experiment with ground station
% 
% ground_LLH = [0.5; -0.3; 1000];
% mainQ2_orbit_determination;
% 
% figure;
% plot(t_24hr_100, initial_error_100);
% grid on;
% hold on;
% plot(t_24hr, initial_error);
% title('Satellite position error: initial orbit determination');
% xlabel('Time (s)');
% ylabel('Position error (m)');
% legend('GS1','GS2', 'Location', 'southeast');
% 
% if save_figures == true
%     saveas(gcf, 'Q2_error_initial_gs.png');
% end
% 
% figure;
% plot(t_24hr_100, refined_error_100);
% grid on;
% hold on;
% plot(t_24hr, refined_error);
% legend('GS1','GS2', 'Location', 'southeast');
% title('Satellite position error: refined orbit determination');
% xlabel('Time (s)');
% ylabel('Position  error (m)');
% 
% if save_figures == true
%     saveas(gcf, 'Q2_error_refined_gs.png');
% end