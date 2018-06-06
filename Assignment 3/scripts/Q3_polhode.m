%% 440305585
% AERO4701
% Assignment 3
%
% Q3_polhode.m

%% define the two ellipsoids

% get the semi-axis lengths for each ellipsoid
energy_ellipsoid = get_energy_ellipsoid(E.total(1), I);
momentum_ellipsoid = get_momentum_ellipsoid(L.total(1), I);

% get limits for plotting
limits = 1.1 * max([energy_ellipsoid.a, energy_ellipsoid.b, energy_ellipsoid.c ...
    momentum_ellipsoid.a, momentum_ellipsoid.b, momentum_ellipsoid.c]);

figure;
p_energy = surf(energy_ellipsoid.x, energy_ellipsoid.y, energy_ellipsoid.z, 'EdgeColor', 'none', 'FaceAlpha', face_alpha, 'FaceColor', 'r');
hold on;
grid on;
axis equal;
xlim([-limits, limits]);
ylim([-limits, limits]);
zlim([-limits, limits]);
p_momentum = surf(momentum_ellipsoid.x, momentum_ellipsoid.y, momentum_ellipsoid.z, 'EdgeColor', 'none', 'FaceAlpha', face_alpha, 'FaceColor', 'b');
p_w = plot3(w.x, w.y, w.z, 'k', 'LineWidth', line_width);
plot3(-w.x, w.y, w.z, 'k', 'LineWidth', line_width);
plot3(w.x, -w.y, w.z, 'k', 'LineWidth', line_width);
plot3(w.x, w.y, -w.z, 'k', 'LineWidth', line_width);


