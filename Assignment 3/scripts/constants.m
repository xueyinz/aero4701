%% 440305585
% AERO4701
% Assignment 3
%
% constants.m

%% rectangular prism parameters

a = 0.23;           % length a [m]
b = 0.01;           % height b [m]
c = 0.11;           % width c [m]
m = 1;              % mass [kg]
w.x = 0.08;         % initial rotation rate [rad/s]
w.y = 0.08;         % initial rotation rate [rad/s]
w.z = pi/2;         % initial rotation rate [rad/s]

%% simulation parameters

t_end = 60;         % total simulation time [s]
dt = 0.001;         % time step [s]
t = (0:dt:t_end)';  % time vector

