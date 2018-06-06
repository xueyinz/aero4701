%% 440305585
% AERO4701
% Assignment 3
%
% Q1_calculations.m

%% initialisations based on the constants that were loaded

t_vector = 0:dt:t_end;          % time vector [s]
num_steps = length(t_vector);   % number of time steps

% initial corner positions and geometry of the rectangular prism (before rotation)
shape.XYZ = [shape.c/2, shape.b/2, shape.a/2];
shape.XYz = [shape.c/2, shape.b/2, -shape.a/2];
shape.XyZ = [shape.c/2, -shape.b/2, shape.a/2];
shape.Xyz = [shape.c/2, -shape.b/2, -shape.a/2];
shape.xYZ = [-shape.c/2, shape.b/2, shape.a/2];
shape.xYz = [-shape.c/2, shape.b/2, -shape.a/2];
shape.xyZ = [-shape.c/2, -shape.b/2, shape.a/2];
shape.xyz = [-shape.c/2, -shape.b/2, -shape.a/2];
shape.faces = [1 2 4 3; 1 3 7 5; 5 6 8 7; 2 4 8 6];
shape.vertices = [shape.XYZ; shape.XYz; shape.XyZ; shape.Xyz; shape.xYZ; shape.xYz; shape.xyZ; shape.xyz];
shape.prism_limits = 0.6*max([shape.a shape.b shape.c]);

% moments of inertia
I.xx = shape.m*(shape.a^2 + shape.b^2)/12;
I.yy = shape.m*(shape.a^2 + shape.c^2)/12;
I.zz = shape.m*(shape.b^2 + shape.c^2)/12;

%% pre-allocate memory for arrays and structures

% w = angular velocity [rad/s]
w.x = NaN(1, num_steps);
w.y = NaN(1, num_steps);
w.z = NaN(1, num_steps);
w = repmat(w, 1, 3);

% wdot = change in angular velocity [rad/s]
wdot.x = NaN(1, num_steps);
wdot.y = NaN(1, num_steps);
wdot.z = NaN(1, num_steps);
wdot = repmat(wdot, 1, 3);

% L = angular momentum [kg.m^2/s]
L.x = NaN(1, num_steps);
L.y = NaN(1, num_steps);
L.z = NaN(1, num_steps);
L.total = NaN(1, num_steps);
L = repmat(L, 1, 3);

% c = constant coefficients for change in angular velocity
c.x = NaN;
c.y = NaN;
c.z = NaN;
c = repmat(c, 1, 3);

% E = rotational kinetic energy
E.xx = NaN(1, num_steps);
E.yy = NaN(1, num_steps);
E.zz = NaN(1, num_steps);
E.total = NaN(1, num_steps);
E = repmat(E, 1, 3);

% absolute body angles
body_angles.psi = NaN(1, num_steps);
body_angles.theta = NaN(1, num_steps);
body_angles.phi = NaN(1, num_steps);
body_angles = repmat(body_angles, 1, 3);

%% loop three times for rotation about the three axes

for ii = 1:3

    % first time step
    t = 1;
    
    % angular momentum
    L(ii).x(t) = I.xx * w_initial.x;
    L(ii).y(t) = I.yy * w_initial.y;
    L(ii).z(t) = I.zz * w_initial.z;
    L(ii).total(t) = norm([L(ii).x(t), L(ii).y(t), L(ii).z(t)]);

    % angular velocity
    w(ii).x(t) = w_initial.x;
    w(ii).y(t) = w_initial.y;
    w(ii).z(t) = w_initial.z;

    % constant coefficients for change in angular velocity
    c(ii).x = (I.yy - I.zz)/I.xx;
    c(ii).y = (I.zz - I.xx)/I.yy;
    c(ii).z = (I.xx - I.yy)/I.zz;

    % angular velocity
    wdot(ii).x(t) = c(ii).x * w(ii).y(t) * w(ii).z(t);
    wdot(ii).y(t) = c(ii).y * w(ii).z(t) * w(ii).x(t);
    wdot(ii).z(t) = c(ii).z * w(ii).x(t) * w(ii).y(t);

    % calculate angular velocity and angular momentum at each time step,
    % skipping the first time step
    for t = 2:num_steps

        % get new angular velocity
        w(ii).x(t) = w(ii).x(t-1) + wdot(ii).x(t-1) * dt;
        w(ii).y(t) = w(ii).y(t-1) + wdot(ii).y(t-1) * dt;
        w(ii).z(t) = w(ii).z(t-1) + wdot(ii).z(t-1) * dt;

        % get new change in angular velocity
        wdot(ii).x(t) = c(ii).x * w(ii).y(t) * w(ii).z(t);
        wdot(ii).y(t) = c(ii).y * w(ii).z(t) * w(ii).x(t);
        wdot(ii).z(t) = c(ii).z * w(ii).x(t) * w(ii).y(t);

        % get new angular momentum
        L(ii).x(t) = I.xx * w(ii).x(t);
        L(ii).y(t) = I.yy * w(ii).y(t);
        L(ii).z(t) = I.zz * w(ii).z(t);
        L(ii).total(t) = norm([L(ii).x(t), L(ii).y(t), L(ii).z(t)]);

    end

    % angles in the inertial frame of reference
    body_angles(ii) = get_body_frame_angles(w(ii).x, w(ii).y, w(ii).z, dt);

    % (d) rotational kinetic energy
    E(ii).xx = (1/2) * I.xx * w(ii).x.^2;
    E(ii).yy = (1/2) * I.yy * w(ii).y.^2;
    E(ii).zz = (1/2) * I.zz * w(ii).z.^2;
    E(ii).total = E(ii).xx + E(ii).yy + E(ii).zz;
    
    % rotate the axes of the initial angular velocities
    w_initial = rotate_initial_angular_velocities(w_initial);

end
