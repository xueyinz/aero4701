%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PlotEarthSphere.m - Script for plotting an Earth sized sphere to use for
% visualisation of you satellite orbit (AERO4701, 2018, Assignment 1). You
% may either use this script or cut and paste the code into your own
% script.
%
% By Mitch Bryson, 2018.
%

% Create figure and load topographical Earth map
figure
load('topo.mat','topo','topomap1');

% Create a sphere, make it earth sized (in meters)
[x,y,z] = sphere(50);
x = x.*6378000;
y = y.*6378000;
z = z.*6378000;

% Set visual properties for plot
props.AmbientStrength = 0.1;
props.DiffuseStrength = 1;
props.SpecularColorReflectance = .5;
props.SpecularExponent = 20;
props.SpecularStrength = 1;
props.FaceColor= 'texture';
props.EdgeColor = 'none';
props.FaceLighting = 'phong';
props.Cdata = topo;

% Plot Earth
surface(x,y,z,props);
title('ECI Co-ordinates of Orbit')
axis equal
grid on
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
hold on
