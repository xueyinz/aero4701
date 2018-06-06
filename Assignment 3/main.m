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

constants_type = '';

while isempty(constants_type)
    
    constants_type = input('Enter ''d'' for default connstants, ''u'' for user constants: ', 's');
    
    
    if isempty(constants_type)
        fprintf('Invalid input\n');
    elseif ((constants_type == 'd') || (constants_type == 'D'))
        fprintf('\nUsing default constants from assignment sheet\n\n');
        default_constants;
        Q1_calculations;
        Q1_answers;
        Q2_animate;
        Q3_polhode;
    elseif ((constants_type == 'u') || (constants_type == 'U'))
        fprintf('\nUsing the user constants defined in <user_constants.m>\n\n');
        user_constants;        
    else
        constants_type = '';
        fprintf('Invalid input\n');
    end
    
end
