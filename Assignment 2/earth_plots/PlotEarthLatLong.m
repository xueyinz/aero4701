%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlotEarthLatLong.m - Script for plotting a Lat-Long map of the Earth to
% use for visualisation of the ground trace of your satellite orbit
% (AERO4701, 2018, Assignment 1). You may either use this script or cut and
% paste the code into your own script.
%
% By Mitch Bryson, 2018.
%

% Create figure and load topographical Earth map
figure
em = imread('earthmap1.png');

% Plot Earth map
image([-180,180],[90,-90],em);
axis xy
% axis equal
title('Ground Trace of Orbit')
grid on
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
hold on
