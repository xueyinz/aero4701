%% 440305585
% AERO4701
% Assignment 3
%
% main.m

clearvars;
close all;
Q1 = true;
Q2 = true;
Q3 = false;

%% initialisation

addpath('./scripts/');

constants;

if Q1
    mainQ1;
%     if Q2
%         rotation_matrix = get_rotation_maxtrix
%     end
end

if Q3
    mainQ3;
end