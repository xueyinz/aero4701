%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1_A_plots.m

PlotEarthSphere;
view(3);

% get all the plot handles for each satellite
for ii = 1:size(orbits,2)
    plot3(orbits(ii).pos_ECI(1,:), orbits(ii).pos_ECI(2,:), orbits(ii).pos_ECI(3,:), 'LineWidth', 0.5, 'LineStyle', '--', 'Color', 'k');
    orbits(ii).plot_handle = plot3(orbits(ii).pos_ECI(1,1), orbits(ii).pos_ECI(2,1), orbits(ii).pos_ECI(3,1), 'o', 'MarkerEdgeColor','g', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
    hold on;
end

% loop through each time step and update position of each satellite
for ii = 1:size(t_12hr,2)
    for jj = 1:size(orbits,2)
        orbits(jj).plot_handle.XData = orbits(jj).pos_ECI(1,ii);
        orbits(jj).plot_handle.YData = orbits(jj).pos_ECI(2,ii);
        orbits(jj).plot_handle.ZData = orbits(jj).pos_ECI(3,ii);
    end
    drawnow;
end