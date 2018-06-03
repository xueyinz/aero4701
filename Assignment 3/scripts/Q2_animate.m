%% 440305585
% AERO4701
% Assignment 3
%
% Q2_animate.m

%% Question 2

% corners of the rectangle
shape.XYZ = NaN(num_steps, 3);
shape.XYZ = NaN(num_steps, 3);
shape.XYz = NaN(num_steps, 3);
shape.XyZ = NaN(num_steps, 3);
shape.Xyz = NaN(num_steps, 3);
shape.xYZ = NaN(num_steps, 3);
shape.xYz = NaN(num_steps, 3);
shape.xyZ = NaN(num_steps, 3);
shape.xyz = NaN(num_steps, 3);
XYZ = [shape.c/2, shape.b/2, shape.a/2];
XYz = [shape.c/2, shape.b/2, -shape.a/2];
XyZ = [shape.c/2, -shape.b/2, shape.a/2];
Xyz = [shape.c/2, -shape.b/2, -shape.a/2];
xYZ = [-shape.c/2, shape.b/2, shape.a/2];
xYz = [-shape.c/2, shape.b/2, -shape.a/2];
xyZ = [-shape.c/2, -shape.b/2, shape.a/2];
xyz = [-shape.c/2, -shape.b/2, -shape.a/2];

% rotation matrix to transform body frame into intertial frame
rotation_matrix = get_rotation_matrix(body_angles.psi, body_angles.theta, body_angles.phi);

% get the corner values at each time step
for t = 1:num_steps
    shape.XYZ(t, :) = XYZ * rotation_matrix(:, :, t);
    shape.XYz(t, :) = XYz * rotation_matrix(:, :, t);
    shape.XyZ(t, :) = XyZ * rotation_matrix(:, :, t);
    shape.Xyz(t, :) = Xyz * rotation_matrix(:, :, t);
    shape.xYZ(t, :) = xYZ * rotation_matrix(:, :, t);
    shape.xYz(t, :) = xYz * rotation_matrix(:, :, t);
    shape.xyZ(t, :) = xyZ * rotation_matrix(:, :, t);
    shape.xyz(t, :) = xyz * rotation_matrix(:, :, t);
end

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

%% animating prism in the body frame

faces = [1 2 4 3; 1 3 7 5; 5 6 8 7; 2 4 8 6];
vertices = [XYZ; XYz; XyZ; Xyz; xYZ; xYz; xyZ; xyz];
limits = 0.6*max([shape.a shape.b shape.c]);
dominant_rotation = get_dominant_rotation_axis(I, [w.x_initial, w.y_initial, w.z_initial]);

figure;
rectangular_prism = patch('Faces', faces, 'Vertices', vertices, 'FaceAlpha', 0);
view(3);
axis equal;
grid on;
xlim([-limits limits]);
ylim([-limits limits]);
zlim([-limits limits]);
% gif('yes.gif','DelayTime',0.05,'LoopCount',1,'frame',gcf);
title_string = sprintf('Dominant rotation about %s (inertial frame)(t = %.2f)', dominant_rotation, 0);
plot_title = title(title_string);

for t = 2:animate_speed:num_steps
    rectangular_prism.Vertices = [shape.XYZ(t,:); shape.XYz(t,:); shape.XyZ(t,:); shape.Xyz(t,:); shape.xYZ(t,:); shape.xYz(t,:); shape.xyZ(t,:); shape.xyz(t,:)];
    title_string = sprintf('Dominant rotation about %s (inertial frame)(t = %.2f)', dominant_rotation, t_vector(t));
    plot_title.String = title_string;
    drawnow;
%     gif;
end
