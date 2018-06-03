%% 440305585
% AERO4701
% Assignment 3
%
% get_rotation_matrix.m

function rotation_matrix = get_rotation_matrix(psi, theta, phi)

    c_psi = cos(psi);
    s_psi = sin(psi);
    c_theta = cos(theta);
    s_theta = sin(theta);
    c_phi = cos(phi);
    s_phi = sin(phi);
    
    rotation_matrix = NaN(3,3,length(psi));
    rotation_matrix(1,1,:) = c_psi.*c_phi - c_theta.*s_phi.*s_psi;
    rotation_matrix(1,2,:) = s_phi.*c_psi + c_theta.*c_phi.*s_psi;
    rotation_matrix(1,3,:) = s_psi.*s_theta;
    
    rotation_matrix(2,1,:) = -s_psi.*c_phi - c_theta.*s_phi.*c_psi;
    rotation_matrix(2,2,:) = -s_phi.*s_psi + c_theta.*c_phi.*c_psi;
    rotation_matrix(2,3,:) = c_psi.*s_theta;
    
    rotation_matrix(3,1,:) = s_theta.*s_phi;
    rotation_matrix(3,2,:) = -s_theta.*c_phi;
    rotation_matrix(3,3,:) = c_theta;
    
end