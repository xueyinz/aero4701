%% 440305585
% AERO4701
% Assignment 3
%
% Q2_animate.m

%% pre-allocate memory for arrays and structs

% corners of the rectangular prism at each time step
shape_sim.XYZ = NaN(num_steps, 3);
shape_sim.XYZ = NaN(num_steps, 3);
shape_sim.XYz = NaN(num_steps, 3);
shape_sim.XyZ = NaN(num_steps, 3);
shape_sim.Xyz = NaN(num_steps, 3);
shape_sim.xYZ = NaN(num_steps, 3);
shape_sim.xYz = NaN(num_steps, 3);
shape_sim.xyZ = NaN(num_steps, 3);
shape_sim.xyz = NaN(num_steps, 3);
shape_sim = repmat(shape_sim, 1, 3);

% dominant rotation axis
dominant_rotation = ["", "", ""];

%% loop three times for rotation about the three axes

for ii = 1:3
    
    % rotation matrix to transform body frame into intertial frame
    rotation_matrix = get_rotation_matrix(body_angles(ii).psi, body_angles(ii).theta, body_angles(ii).phi);
    
    % get the corner values at each time step
    for t = 1:num_steps
        shape_sim(ii).XYZ(t, :) = shape.XYZ * rotation_matrix(:, :, t);
        shape_sim(ii).XYz(t, :) = shape.XYz * rotation_matrix(:, :, t);
        shape_sim(ii).XyZ(t, :) = shape.XyZ * rotation_matrix(:, :, t);
        shape_sim(ii).Xyz(t, :) = shape.Xyz * rotation_matrix(:, :, t);
        shape_sim(ii).xYZ(t, :) = shape.xYZ * rotation_matrix(:, :, t);
        shape_sim(ii).xYz(t, :) = shape.xYz * rotation_matrix(:, :, t);
        shape_sim(ii).xyZ(t, :) = shape.xyZ * rotation_matrix(:, :, t);
        shape_sim(ii).xyz(t, :) = shape.xyz * rotation_matrix(:, :, t);
    end
    
    % animating prism in the body frame
    dominant_rotation(ii) = get_dominant_rotation_axis(I, [w_initial.x, w_initial.y, w_initial.z]);
    
    % rotate the axes of the initial angular velocities
    w_initial = rotate_initial_angular_velocities(w_initial);

end
