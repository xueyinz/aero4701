%% 440305585
% AERO4701
% Assignment 3
%
% plot_everything.m

%% plotting angular velocity

figure;
plot(t_vector, w.x);
hold on;
grid on;
plot(t_vector, w.y);
plot(t_vector, w.z);
title('Angular velocity');

%% plotting angular momentum

figure;
plot(t_vector, L.x);
hold on;
grid on;
plot(t_vector, L.y);
plot(t_vector, L.z);
plot(t_vector, L.total);
title('Angular momentum');

%% plotting rotational kinetic energy

figure;
plot(t_vector, E.xx);
hold on;
grid on;
plot(t_vector, E.yy);
plot(t_vector, E.zz);
plot(t_vector, E.total);
title('Rotational kinetic energy');

%% prism rotating in the inertial frame

figure;
rectangular_prism = patch('Faces', faces, 'Vertices', vertices, 'FaceAlpha', 0);
view(3);
axis equal;
grid on;
xlim([-prism_limits prism_limits]);
ylim([-prism_limits prism_limits]);
zlim([-prism_limits prism_limits]);
title_string = sprintf('Dominant rotation about %s (inertial frame)(t = %.2f)', dominant_rotation, 0);
plot_title = title(title_string);

for t = 2:animate_speed:num_steps
    rectangular_prism.Vertices = [shape.XYZ(t,:); shape.XYz(t,:); shape.XyZ(t,:); shape.Xyz(t,:); shape.xYZ(t,:); shape.xYz(t,:); shape.xyZ(t,:); shape.xyz(t,:)];
    title_string = sprintf('Dominant rotation about %s (inertial frame)(t = %.2f)', dominant_rotation, t_vector(t));
    plot_title.String = title_string;
    drawnow;
end

%% Poinsot polhode

figure;
p_energy = surf(energy_ellipsoid.x, energy_ellipsoid.y, energy_ellipsoid.z, 'EdgeColor', 'none', 'FaceAlpha', face_alpha, 'FaceColor', 'r');
hold on;
grid on;
axis equal;
xlim([-ellipsoid_limits, ellipsoid_limits]);
ylim([-ellipsoid_limits, ellipsoid_limits]);
zlim([-ellipsoid_limits, ellipsoid_limits]);
p_momentum = surf(momentum_ellipsoid.x, momentum_ellipsoid.y, momentum_ellipsoid.z, 'EdgeColor', 'none', 'FaceAlpha', face_alpha, 'FaceColor', 'b');
p_w = plot3(w.x, w.y, w.z, 'k', 'LineWidth', line_width);
plot3(-w.x, w.y, w.z, 'k', 'LineWidth', line_width);
plot3(w.x, -w.y, w.z, 'k', 'LineWidth', line_width);
plot3(w.x, w.y, -w.z, 'k', 'LineWidth', line_width);
p_w_point = scatter3(w.x(1), w.y(2), w.z(3), 20*line_width, 'y', 'filled');
for t = 2:animate_speed:num_steps
    p_w_point.XData = w.x(t);
    p_w_point.YData = w.y(t);
    p_w_point.ZData = w.z(t);
    drawnow;
end



