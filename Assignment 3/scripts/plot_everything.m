%% 440305585
% AERO4701
% Assignment 3
%
% plot_everything.m

plot_id = 1;

%% plotting angular velocity

figure(plot_id);
plot_id = plot_id + 1;
for ii = 1:3
    
    subplot(3, 1, ii);
    
    plot(t_vector, w(ii).x, 'LineWidth', line_width);
    hold on;
    grid on;
    plot(t_vector, w(ii).y, 'LineWidth', line_width);
    plot(t_vector, w(ii).z, 'LineWidth', line_width);
    
    title_string = sprintf('Main rotation about %s (%s-axis)',...
        dominant_rotation(ii).I, dominant_rotation(ii).axis);
    title(title_string);
    legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'eastoutside');

end
saveas(gcf, 'angular_velocity.png');

%% plotting angular momentum

figure(plot_id);
plot_id = plot_id + 1;
for ii = 1:3
    
    subplot(3, 1, ii);

    plot(t_vector, L(ii).x, 'LineWidth', line_width);
    hold on;
    grid on;
    plot(t_vector, L(ii).y, 'LineWidth', line_width);
    plot(t_vector, L(ii).z, 'LineWidth', line_width);
    plot(t_vector, L(ii).total, 'LineWidth', line_width);
    
    title_string = sprintf('Main rotation about %s (%s-axis)',...
        dominant_rotation(ii).I, dominant_rotation(ii).axis);
    title(title_string);
    legend('L_x', 'L_y', 'L_z', 'L_total', 'Location', 'eastoutside');

end
saveas(gcf, 'angular_momentum.png');

%% plotting rotational kinetic energy

figure(plot_id);
plot_id = plot_id + 1;
for ii = 1:3
    
    subplot(3, 1, ii);

    plot(t_vector, E(ii).xx, 'LineWidth', line_width);
    hold on;
    grid on;
    plot(t_vector, E(ii).yy, 'LineWidth', line_width);
    plot(t_vector, E(ii).zz, 'LineWidth', line_width);
    plot(t_vector, E(ii).total, 'LineWidth', line_width);

    title_string = sprintf('Main rotation about %s (%s-axis)',...
        dominant_rotation(ii).I, dominant_rotation(ii).axis);
    title(title_string);
    legend('E_x', 'E_y', 'E_z', 'E_total', 'Location', 'eastoutside');
    
end
saveas(gcf, 'kinetic_energy.png');

%% animations

% colour of the edges in the prism
col = [0,     0.4470, 0.7410;
    0.8500, 0.3250, 0.0980;
    0.9290, 0.6940, 0.1250;
    0.4940, 0.1840, 0.5560];

% create the plots for the first time step
for ii = 1:3
    
    figure(plot_id);
    plot_id = plot_id + 1;
    
    % prism in inertial frame
    subplot(1, 2, 1);
    rectangular_prism = patch('Faces', shape.faces, 'Vertices', shape.vertices,...
        'FaceVertexCData', col, 'FaceColor', 'flat');
    view(3);
    grid on;
    hold on;
    axis equal;
    xlim([-shape.prism_limits, shape.prism_limits]);
    ylim([-shape.prism_limits, shape.prism_limits]);
    zlim([-shape.prism_limits, shape.prism_limits]);
    title('Inertial frame');
    scatter3(0, 0, 0, 20*line_width, 'k', 'filled');
    
    % Poinsot polhode
    subplot(1, 2, 2);
    p_energy = surf(energy_ellipsoid(ii).x, energy_ellipsoid(ii).y, energy_ellipsoid(ii).z,...
        'EdgeColor', 'none', 'FaceAlpha', face_alpha, 'FaceColor', 'r');
    grid on;
    hold on;
    axis equal;
    xlim([-ellipsoid_limits(ii), ellipsoid_limits(ii)]);
    ylim([-ellipsoid_limits(ii), ellipsoid_limits(ii)]);
    zlim([-ellipsoid_limits(ii), ellipsoid_limits(ii)]);
    title('Omega space (polhode)');
    p_momentum = surf(momentum_ellipsoid(ii).x, momentum_ellipsoid(ii).y, momentum_ellipsoid(ii).z,...
        'EdgeColor', 'none', 'FaceAlpha', face_alpha, 'FaceColor', 'b');
    p_w = plot3(w(ii).x, w(ii).y, w(ii).z, 'y', 'LineWidth', line_width);
    plot3(-w(ii).x, w(ii).y, w(ii).z, 'y', 'LineWidth', line_width);
    plot3(w(ii).x, -w(ii).y, w(ii).z, 'y', 'LineWidth', line_width);
    plot3(w(ii).x, w(ii).y, -w(ii).z, 'y', 'LineWidth', line_width);
    p_w_point = scatter3(w(ii).x(1), w(ii).y(1), w(ii).z(1), 20*line_width, 'k', 'filled');
    p_legend = legend([p_energy, p_momentum, p_w, p_w_point],...
        'Energy ellipsoid', 'Angular momentum ellipsoid', 'Intersection of ellipsoids', 'Current angular momentum vector',...
        'Location', 'south');
    set(p_legend, 'Position', [0.572857142857143 0.0933333333333334 0.321428571428571 0.117857142857143]);
    
    % title for whole figure
    title_string = sprintf('Main rotation about %s (%s-axis) (t = 0s)',...
        dominant_rotation(ii).I, dominant_rotation(ii).axis);
    a = annotation('textbox', [0 0.8 1 0.1], ...
    'String', title_string, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'FontSize', 14);
    
    % update the frame for each time step
    for t = 2:animate_speed:num_steps
        
        % title
        a.String = sprintf('Main rotation about %s (%s-axis) (t = %.1fs)',...
            dominant_rotation(ii).I, dominant_rotation(ii).axis, t_vector(t));

        % update angular velocity vector
        subplot(1, 2, 1);
        p_w_point.XData = w(ii).x(t);
        p_w_point.YData = w(ii).y(t);
        p_w_point.ZData = w(ii).z(t);

        % update rectangular prism
        subplot(1, 2, 2);
        rectangular_prism.Vertices = [shape_sim(ii).XYZ(t,:); shape_sim(ii).XYz(t,:);...
            shape_sim(ii).XyZ(t,:); shape_sim(ii).Xyz(t,:); shape_sim(ii).xYZ(t,:);...
            shape_sim(ii).xYz(t,:); shape_sim(ii).xyZ(t,:); shape_sim(ii).xyz(t,:)];

        drawnow;

    end

end

