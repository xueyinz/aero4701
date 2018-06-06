%% 440305585
% AERO4701
% Assignment 3
%
% rotate_initial_angular_velocities.m

function w_initial_out = rotate_initial_angular_velocities(w_initial)

    w_initial_out.x = w_initial.z;
    w_initial_out.y = w_initial.x;
    w_initial_out.z = w_initial.y;
    
end
