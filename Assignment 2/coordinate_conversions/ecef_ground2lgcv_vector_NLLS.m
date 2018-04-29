%% 440305585
% AERO4701
%
% Convert Earth-Centred Earth-Fixed coordinates of a satellite to local 
% geocentric vertical coordinates relative to a ground station in latitude,
% longitude, and height.
%
% inputs:   pos_satellite_ECEF = [x; y; z] = ECEF coordinates [m, m, m] (x, y, z are vectors)
%           pos_ground_ECEF = [x; y; z] = ground station in ECEF coordinates [m, m, m]
% outputs:  pos_LGCV = [x; y; z] = LGCV coordinates [m, m, m]

function pos_satellite_relative_LGCV = ecef_ground2lgcv_vector_NLLS(pos_satellite_ECEF, pos_ground_LLH, pos_ground_ECEF)
    
    % convert groundstation geocentric LLH to ECEF
    pos_relative = pos_satellite_ECEF - pos_ground_ECEF;
    
    lambda = pos_ground_LLH(1);      % latitude
    phi = pos_ground_LLH(2);         % longitude
    
    s_lambda = sin(lambda);
    c_lambda = cos(lambda);
    s_phi = sin(phi);
    c_phi = cos(phi);
    
    % rotation matrix
    C = [-s_lambda*c_phi,   -s_phi,     -c_lambda*c_phi;
        -s_lambda*s_phi,    c_phi,      -c_lambda*s_phi;
        c_lambda,           0,          -s_lambda];
    
    pos_satellite_relative_LGCV = C\pos_relative;

end