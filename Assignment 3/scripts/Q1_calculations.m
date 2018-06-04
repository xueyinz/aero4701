%% 440305585
% AERO4701
% Assignment 3
%
% Q1_calculations.m

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

if question == 1
    
    E.xx = (1/2) * I.xx * w.x.^2;
    E.yy = (1/2) * I.yy * w.y.^2;
    E.zz = (1/2) * I.zz * w.z.^2;
    E.total = E.xx + E.yy + E.zz;

    %% (e) copper wire temperature

    % angular momentum is conserved for the system
    w.final = L.total_initial/max([I.xx I.yy I.zz]);
    
    % rotational kinetic energy in the final stable state
    E.minimum = (1/2) * (max([I.xx I.yy I.zz]) * w.final^2);
    % E.minimum = (1/2) * L.total_initial^2 / min([I.xx, I.yy, I.zz]);
    
    % loss in rotational kinetic energy
    E.loss = E.total(1) - E.minimum;
    
    % wire properties
    wire.volume = pi * (wire.d/2)^2 * wire.l;
    wire.mass = wire.density * wire.volume;             % [g]
    
    % change in thermal energy for the wire
    wire.delta_temp = E.loss / (wire.mass * wire.c);
%     wire.final_temp = wire.temp + wire.delta_temp;

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
    fprintf('\tE_rotational = %f J\n', E.total(1));

    fprintf('1. (e)\n');
    fprintf('\tChange in temperature of wire = %f celsius\n', wire.delta_temp);
    
end
