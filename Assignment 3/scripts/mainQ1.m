%% 440305585
% AERO4701
% Assignment 3
%
% mainQ1.m

%% initialisations

t_vector = 0:dt:t_end;          % time vector [s]
num_steps = length(t_vector);   % number of time steps

% faster to operate on a pre-allocated array
w.x = NaN(1, num_steps);        % x angular velocity vector [rad/s]
w.y = NaN(1, num_steps);        % y angular velocity vector [rad/s]
w.z = NaN(1, num_steps);        % z angular velocity vector [rad/s]

wdot.x = NaN(1, num_steps);     % x change in angular velocity vector [rad/s]
wdot.y = NaN(1, num_steps);     % y change in angular velocity vector [rad/s]
wdot.z = NaN(1, num_steps);     % z change in angular velocity vector [rad/s]

L.x = NaN(1, num_steps);        % x angular momentum vector [kg.m^2/s]
L.y = NaN(1, num_steps);        % y angular momentum vector [kg.m^2/s]
L.z = NaN(1, num_steps);        % z angular momentum vector [kg.m^2/s]
L.total = NaN(1, num_steps);    % total angular momentum vector [kg.m^2/s]



%% (a) moments of inertia

I.xx = shape.m*(shape.a^2 + shape.b^2)/12;
I.yy = shape.m*(shape.a^2 + shape.c^2)/12;
I.zz = shape.m*(shape.b^2 + shape.c^2)/12;

%% (b) inertial angular momentum magnitude

L.x_initial = I.xx * w.x_initial;
L.y_initial = I.yy * w.y_initial;
L.z_initial = I.zz * w.z_initial;
L.total_initial = norm([L.x_initial, L.y_initial, L.z_initial]);

%% (c) maximum angular velocity

% initialise the first time step
t = 1;

% angular velocity
w.x(t) = w.x_initial;
w.y(t) = w.y_initial;
w.z(t) = w.z_initial;

% constant coefficients for change in angular velocity
c.x = (I.yy - I.zz)/I.xx;
c.y = (I.zz - I.xx)/I.yy;
c.z = (I.xx - I.yy)/I.zz;

% angular velocity
wdot.x(t) = c.x * w.y(t) * w.z(t);
wdot.y(t) = c.y * w.z(t) * w.x(t);
wdot.z(t) = c.z * w.x(t) * w.y(t);

% angular momentum
L.x(t) = L.x_initial;
L.y(t) = L.y_initial;
L.z(t) = L.z_initial;
L.total(t) = L.total_initial;

% calculate angular velocity and angular momentum at each time step,
% skipping the first time step
for t = 2:num_steps
    
    % get new angular velocity
    w.x(t) = w.x(t-1) + wdot.x(t-1) * dt;
    w.y(t) = w.y(t-1) + wdot.y(t-1) * dt;
    w.z(t) = w.z(t-1) + wdot.z(t-1) * dt;
    
    % get new change in angular velocity
    wdot.x(t) = c.x * w.y(t) * w.z(t);
    wdot.y(t) = c.y * w.z(t) * w.x(t);
    wdot.z(t) = c.z * w.x(t) * w.y(t);
    
    % get new angular momentum
    L.x(t) = I.xx * w.x(t);
    L.y(t) = I.yy * w.y(t);
    L.z(t) = I.zz * w.z(t);
    L.total(t) = norm([L.x(t), L.y(t), L.z(t)]);
    
end

% maximum angular velocity
w.x_max = max(w.x);
w.y_max = max(w.y);
w.z_max = max(w.z);

% angles in the inertial frame of reference
body_angles = get_body_frame_angles(w.x, w.y, w.z, dt);

%% (d) rotational kinetic energy

E.xx = (1/2) * I.xx * w.x.^2;
E.yy = (1/2) * I.yy * w.y.^2;
E.zz = (1/2) * I.zz * w.z.^2;
E.total = E.xx + E.yy + E.zz;

%% (e) copper wire temperature

% E.minimum = (1/2) * L.total_initial^2 / min([I.xx, I.yy, I.zz]);
% E.diff = E.total(1) - E.minimum;

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

%% results

fprintf('1. (a)\n');
fprintf('\tI_xx = %f kg.m^2\n', I.xx);
fprintf('\tI_yy = %f kg.m^2\n', I.yy);
fprintf('\tI_zz = %f kg.m^2\n', I.zz);

fprintf('1. (b)\n');
fprintf('\tL_total = %f kg.m^2/s\n', L.total_initial);

fprintf('1. (c)\n');
fprintf('\tmax(w_x) = %f rad/s\n', w.x_max);
fprintf('\tmax(w_y) = %f rad/s\n', w.y_max);
fprintf('\tmax(w_z) = %f rad/s\n', w.z_max);

fprintf('1. (d)\n');
fprintf('\tmax(w_x) = %f rad/s\n', w.x_max);

%% plotting

% angular velocity
figure;
plot(t_vector, w.x);
hold on;
grid on;
plot(t_vector, w.y);
plot(t_vector, w.z);
title('Angular velocity');

% angular momentum
figure;
plot(t_vector, L.x);
hold on;
grid on;
plot(t_vector, L.y);
plot(t_vector, L.z);
plot(t_vector, L.total);
title('Angular momentum');

% energy
figure;
plot(t_vector, E.xx);
hold on;
grid on;
plot(t_vector, E.yy);
plot(t_vector, E.zz);
plot(t_vector, E.total);
title('Energy');

% body frame
faces = [1 2 4 3; 1 3 7 5; 5 6 8 7; 2 4 8 6];
vertices = [XYZ; XYz; XyZ; Xyz; xYZ; xYz; xyZ; xyz];
limits = 0.6*max([shape.a shape.b shape.c]);

figure;
rectangular_prism = patch('Faces', faces, 'Vertices', vertices, 'FaceAlpha', 0);
view(3);
axis equal;
grid on;
xlim([-limits limits]);
ylim([-limits limits]);
zlim([-limits limits]);

for t = 2:animate_speed:num_steps
    rectangular_prism.Vertices = [shape.XYZ(t,:); shape.XYz(t,:); shape.XyZ(t,:); shape.Xyz(t,:); shape.xYZ(t,:); shape.xYz(t,:); shape.xyZ(t,:); shape.xyz(t,:)];
    drawnow;
end

