%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1_A_plots.m

PlotEarthSphere;
view(3);

% plot all the orbits
for ii = 1:n_satellites
    plot3(orbits(ii).pos_ECI(1,:), orbits(ii).pos_ECI(2,:), orbits(ii).pos_ECI(3,:));
    hold on;
end

if save_figures == true
    saveas(gcf, 'Q1A_12hr_orbits.png');
end

% get all the plot handles for each satellite
for ii = 1:n_satellites
    orbits(ii).plot_handle = plot3(orbits(ii).pos_ECI(1,1), orbits(ii).pos_ECI(2,1), orbits(ii).pos_ECI(3,1), 'ko');
end

% loop through each time step and update position of each satellite
for ii = 1:size(t_12hr,2)
    for jj = 1:n_satellites
        orbits(jj).plot_handle.XData = orbits(jj).pos_ECI(1,ii);
        orbits(jj).plot_handle.YData = orbits(jj).pos_ECI(2,ii);
        orbits(jj).plot_handle.ZData = orbits(jj).pos_ECI(3,ii);
    end
    drawnow;
end

if save_figures == true
    saveas(gcf, 'Q1A_12hr_orbits_instantaneous.png');
end