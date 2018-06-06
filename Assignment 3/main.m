%% 440305585
% AERO4701
% Assignment 3
%
% main.m

clearvars;
close all;

%% initialisation

addpath('./scripts/', './functions/');

default_constants;

%% Question 1

question = 1;
Q1_calculations;

%% Question 2

question = 2;
Q1_calculations;
% Q2_animate;
% 
% % y-axis
% w.x_initial = 0.08;             % initial x angular velocity [rad/s]
% w.y_initial = pi/2;             % initial y angular velocity [rad/s]
% w.z_initial = 0.1;              % initial z angular velocity [rad/s]
% Q1_calculations;
% Q2_animate;
% 
% % z-axis
% w.x_initial = pi/2;             % initial x angular velocity [rad/s]
% w.y_initial = 0.06;             % initial y angular velocity [rad/s]
% w.z_initial = 0.08;             % initial z angular velocity [rad/s]
% Q1_calculations;
% Q2_animate;

%% Question 3

Q3_polhode;
