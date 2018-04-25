%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1_B_plots.m

figure;
plot3(uav_truth.LGCV(1,:), uav_truth.LGCV(2,:), uav_truth.LGCV(3,:),'.');
hold on;
plot3(uav.pos_LGCV(1,:), uav.pos_LGCV(2,:), uav.pos_LGCV(3,:), '.');

title('UAV in LGCV coordinates');
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
grid on;