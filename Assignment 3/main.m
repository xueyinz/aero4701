%% 440305585
% AERO4701
% Assignment 3
%
% main.m
%
% Note to marker: to simulate using different initial rotational
% velocities, slab dimensions, or slab mass, you can modify the code inside
% <user_constants.m>. Enter 'u' at run-time when prompted.

clearvars;
close all;

%% initialisation

addpath('./scripts/', './functions/');

% determines which constants file to load
constants_type = 'd';
        
% check user input
if constants_type == 'd'
    fprintf('\nLoading default constants from assignment sheet...\n\n');
    default_constants;
else
    fprintf('\nLoading constants defined in <user_constants.m>...\n\n');
    user_constants;        
end
    
%% run simulation scripts

% analyse rotational states
Q1_calculations;

% print answers if using default assignment constants
if constants_type == 'd'
    Q1_answers;
end

% prepare animation of rectangular prism in inertial frame
Q2_animate;

% prepare polhode data
Q3_polhode;

% let the simulation run!
saving = false;
saving_type = '';
plot_everything;

% % graphs for report purposes
% plot_report;