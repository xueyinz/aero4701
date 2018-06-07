% for report

close all;

% min(body_angles(1).psi)
% min(body_angles(1).theta)
% min(body_angles(1).phi)
% max(body_angles(1).psi)
% max(body_angles(1).theta)
% max(body_angles(1).phi)
% 
% min(body_angles(2).psi)
% min(body_angles(2).theta)
% min(body_angles(2).phi)
% max(body_angles(2).psi)
% max(body_angles(2).theta)
% max(body_angles(2).phi)
% 
% min(body_angles(3).psi)
% min(body_angles(3).theta)
% min(body_angles(3).phi)
% max(body_angles(3).psi)
% max(body_angles(3).theta)
% max(body_angles(3).phi)

% figure;
% plot(t_vector, w(ii).x, 'LineWidth', line_width);
% hold on;
% grid on;
% plot(t_vector, w(ii).y, 'LineWidth', line_width);
% plot(t_vector, w(ii).z, 'LineWidth', line_width);
% title('Angular velocity \omega vs time');
% legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'eastoutside');
% xlabel('Time (s)');
% ylabel('\omega (rad/s)');
% saveas(gcf, 'angular_velocity.png');
% 
% figure;
% plot(t_vector, L(ii).x, 'LineWidth', line_width);
% hold on;
% grid on;
% plot(t_vector, L(ii).y, 'LineWidth', line_width);
% plot(t_vector, L(ii).z, 'LineWidth', line_width);
% plot(t_vector, L(ii).total, 'LineWidth', line_width);
% title('Angular momentum L vs time');
% legend('L_x', 'L_y', 'L_z', 'L_t_o_t_a_l', 'Location', 'eastoutside');
% xlabel('Time (s)');
% ylabel('L (kg.m^2/s)');
% saveas(gcf, 'angular_momentum.png');
% 
% figure;
% plot(t_vector, E(ii).xx, 'LineWidth', line_width);
% hold on;
% grid on;
% plot(t_vector, E(ii).yy, 'LineWidth', line_width);
% plot(t_vector, E(ii).zz, 'LineWidth', line_width);
% plot(t_vector, E(ii).total, 'LineWidth', line_width);
% title('Rotational kinetic energy E vs time');
% legend('E_x', 'E_y', 'E_z', 'E_t_o_t_a_l', 'Location', 'eastoutside');
% xlabel('Time (s)');
% ylabel('E (J)');
% saveas(gcf, 'kinetic_energy.png');

% create the plots for the first time step
for ii = 1:3
    
    figure(ii);
    
    % prism in inertial frame
    rectangular_prism = patch('Faces', shape.faces, 'Vertices', shape.vertices,...
        'FaceColor', 'w', 'FaceAlpha', 0, 'LineStyle', '--');
    view(3);
    grid on;
    hold on;
    axis equal;
    xlim([-shape.prism_limits, shape.prism_limits]);
    ylim([-shape.prism_limits, shape.prism_limits]);
    zlim([-shape.prism_limits, shape.prism_limits]);
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('z (m)');
    title_string = sprintf('Inertial frame for main rotation about %s (%s-axis)',...
        dominant_rotation(ii).I, dominant_rotation(ii).axis);
    title(title_string);
    scatter3(0, 0, 0, 5*line_width, 'k', 'filled');
    s_corner = scatter3(shape_sim(ii).XYZ(1,1), shape_sim(ii).XYZ(1,2), shape_sim(ii).XYZ(1,3), 20*line_width, 'b', 'filled');
    p_corner = plot3(shape_sim(ii).XYZ(1,1), shape_sim(ii).XYZ(1,2), shape_sim(ii).XYZ(1,3), 'r', 'MarkerSize', 2);
    legend([rectangular_prism, s_corner, p_corner], 'Rectangular slab', 'Instantaneous corner position', 'Corner trajectory', 'Location', 'eastoutside');
    
    figure(ii + 3);
    plot(t_vector, body_angles(ii).psi, 'LineWidth', line_width);
    hold on;
    grid on;
    plot(t_vector, body_angles(ii).theta, 'LineWidth', line_width);
    plot(t_vector, body_angles(ii).phi, 'LineWidth', line_width);
    title_string = sprintf('Body angles for main rotation about %s (%s-axis)',...
        dominant_rotation(ii).I, dominant_rotation(ii).axis);
    title(title_string);
    legend('\psi', '\theta', '\phi', 'Location', 'northwest');
    xlabel('Time (s)');
    ylabel('Body angles (rad)');
    saving_string = sprintf('body_angles_%d.png', ii);
    saveas(gcf, saving_string);
    
    % update the frame for each time step
    for t = 2:animate_speed:num_steps

        figure(ii);
        
        % update rectangular prism
        rectangular_prism.Vertices = [shape_sim(ii).XYZ(t,:); shape_sim(ii).XYz(t,:);...
            shape_sim(ii).XyZ(t,:); shape_sim(ii).Xyz(t,:); shape_sim(ii).xYZ(t,:);...
            shape_sim(ii).xYz(t,:); shape_sim(ii).xyZ(t,:); shape_sim(ii).xyz(t,:)];
        s_corner.XData = shape_sim(ii).XYZ(t,1);
        s_corner.YData = shape_sim(ii).XYZ(t,2);
        s_corner.ZData = shape_sim(ii).XYZ(t,3);
        p_corner.XData = shape_sim(ii).XYZ(1:t,1);
        p_corner.YData = shape_sim(ii).XYZ(1:t,2);
        p_corner.ZData = shape_sim(ii).XYZ(1:t,3);
        
        drawnow;

    end
    
    saving_string = sprintf('inertial_%d.png', ii);
    saveas(gcf, saving_string);

end