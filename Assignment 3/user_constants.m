%% 440305585
% AERO4701
% Assignment 3
%
% user_constants.m
%
% Constants that the user can input.

%% simulation parameters

t_end = 60;                     % total simulation time [s]
dt = 0.001;                     % time step [s]
animate_speed = 100;            % speed multiplier for animation

%% rectangular prism parameters

shape.a = 0.23;                 % length a [m]
shape.b = 0.01;                 % height b [m]
shape.c = 0.11;                 % width c [m]

shape.m = 10;                    % mass [kg]

w_initial.x = 0.06;             % initial x angular velocity [rad/s]
w_initial.y = 0.08;             % initial y angular velocity [rad/s]
w_initial.z = pi/2;             % initial z angular velocity [rad/s]

%% flexible copper wire parameters

wire.l = 17.7;                  % length [cm]
wire.d = 0.1;                   % diameter [cm]
wire.temp = 0;                  % temperature [celsius]
wire.c = 0.385;                 % specific heat [J/g]
wire.density = 9;               % density [g/cm^3]

%% simulation parameters

global surf_size;
surf_size = 49;                 % surf mesh size for the ellipsoids

face_alpha = 0.4;               % transparency for the ellipsoids
line_width = 1.5;               % angular velocity vector line widths
