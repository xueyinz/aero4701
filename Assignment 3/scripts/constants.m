%% 440305585
% AERO4701
% Assignment 3
%
% constants.m

%% simulation parameters

t_end = 60;                     % total simulation time [s]
dt = 0.001;                     % time step [s]
t_vector = 0:dt:t_end;          % time vector [s]
num_steps = length(t_vector);   % number of time steps

%% rectangular prism parameters

a = 0.23;                       % length a [m]
b = 0.01;                       % height b [m]
c = 0.11;                       % width c [m]
m = 1;                          % mass [kg]
w.x_initial = 0.08;             % initial x angular velocity [rad/s]
w.y_initial = 0.08;             % initial y angular velocity [rad/s]
w.z_initial = pi/2;             % initial z angular velocity [rad/s]

%% initialise results

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

%% flexible copper wire parameters


