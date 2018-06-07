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

% ask for user input to determine which constants file to load
constants_type = '';
while isempty(constants_type)
    
    constants_type = input('Enter ''d'' to load default constants, or ''u'' to load user-defined constants: ', 's');
    
    % check user input
    if isempty(constants_type)
        fprintf('Invalid input\n');
    elseif ((constants_type == 'd') || (constants_type == 'D'))
        constants_type = 'd';
        fprintf('\nLoading default constants from assignment sheet...\n\n');
        default_constants;
    elseif ((constants_type == 'u') || (constants_type == 'U'))
        fprintf('\nLoading constants defined in <user_constants.m>...\n\n');
        user_constants;        
    else
        constants_type = '';
        fprintf('Invalid input\n');
    end
    
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
saving = true;
saving_type = 'large_dt';
plot_everything;
