%% 440305585
% AERO4701
% Assignment 2
%
% mainQ2.m

clearvars;
close all;

plotting = true;
save_figures = true;   % used for when saving graphs for the report
report = true;

%% initialisation

addpath('./Q2_scripts/', './coordinate_conversions/', './orbit_simulation/');

% load constants
Q2_constants;

mainQ2_orbit_determination;

%% plot results

if plotting == true
    
    mainQ2_plots;
      
end

%% secondary results for the report

if report == true
    
    mainQ2_secondary;
    
end