%% 440305585
% AERO4701
% Assignment 3
%
% get_angles_body_frame.m
%
% Area under the angular velocity curve gives angle. Use a trapezoidal
% approximation for area under the curve

function body_angles = get_body_frame_angles(w_x, w_y, w_z, dt)

    body_angles.psi = dt*cumtrapz(w_x);
    body_angles.theta = dt*cumtrapz(w_y);
    body_angles.phi = dt*cumtrapz(w_z);
    
end
