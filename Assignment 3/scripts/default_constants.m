%% 440305585
% AERO4701
% Assignment 3
%
% default_constants.m
%
% Default constants specified by the assignment.

%% simulation parameters

t_end = 60;                     % total simulation time [s]
dt = 0.001;                     % time step [s]
animate_speed = 50;             % speed multiplier for animation

%% rectangular prism parameters

shape.a = 0.23;                 % length a [m]
shape.b = 0.01;                 % height b [m]
shape.c = 0.11;                 % width c [m]
shape.m = 1;                    % mass [kg]

w.x_initial = 0.08;             % initial x angular velocity [rad/s]
w.y_initial = 0.06;             % initial y angular velocity [rad/s]
w.z_initial = pi/2;             % initial z angular velocity [rad/s]

%% flexible copper wire parameters

wire.l = 17.7;                  % length [cm]
wire.d = 0.1;                   % diameter [cm]
wire.temp = 0;                  % temperature [celsius]


