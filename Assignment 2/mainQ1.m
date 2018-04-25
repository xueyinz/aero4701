%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1.m

clearvars;
close all;

plotting = true;
save_figures = false;   % used for when saving graphs for the report

%% initialisation

addpath('./Q1_data/', './Q1_useful_code/', './Q2_useful_code/', './coordinate_conversions/', './orbit_simulation/', './earth_plots/');
constants;              % load constants

% load GPS ephemeris data
GPS_ephemeris = load_GPS_ephemeris('GPSsat_ephem.txt');
n_satellites = size(GPS_ephemeris, 1);

% load GPS pseudorange data for UAV
GPS_pseudorange = load_GPS_pseudorange('GPS_pseudorange.txt');  % [ time (s), sat #, pseudorange (s) ]
t_uav = [GPS_pseudorange.time];                                 % vector of unique time-stamps in pseudorange data

% load actual UAV positions for comparison
uav_truth = load_UAV_position('UAVPosition.txt');

%% Q1.A - satellite trajectory over a 12hr period

mainQ1_A;

%% Q1.B - satellite trajectory over the time period of the GPS pseudorange data

mainQ1_B;

%% plot results

if plotting == true
    
%     mainQ1_A_plots;
    
    mainQ1_B_plots;
    
end