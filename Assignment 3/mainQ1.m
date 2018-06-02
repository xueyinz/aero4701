%% 440305585
% AERO4701
% Assignment 3
%
% mainQ1.m

clearvars;
close all;

%% initialisation

addpath('./scripts/');

% load constants
constants;

%% (a) moments of inertia

I.xx = m*(a^2 + b^2)/12;
I.yy = m*(a^2 + c^2)/12;
I.zz = m*(b^2 + c^2)/12;
% I = [I_xx 0 0; 0 I_yy 0; 0 0 I_zz];

%% (b) inertial angular momentum vector

L.x = I.xx*w.x;
L.y = I.yy*w.y;
L.z = I.zz*w.z;
L.total = norm([L.x L.y L.z]);

%% 

%% results

fprintf('1. (a)\n');
fprintf('\tI_xx = %f\n', I.xx);
fprintf('\tI_yy = %f\n', I.yy);
fprintf('\tI_zz = %f\n', I.zz);

fprintf('1. (b)\n');
fprintf('\tL_total = %f\n', L.total);
