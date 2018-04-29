%% 440305585
% AERO4701
%
% Gets transformation matrix to convert LGCV into ECEF coordinates.
%
% inputs:   ground_LLH = [x; y; z] = ground station in LLH coordinates [m, m, m]
% outputs:  C = 3x3 transformation matrix

function C = C_lgcv2ecef(ground_LLH)

    lambda = ground_LLH(1);      % latitude
    phi = ground_LLH(2);         % longitude
    
    % rotation matrix
    C = [-sin(lambda)*cos(phi), -sin(phi),  -cos(lambda)*cos(phi);
        -sin(lambda)*sin(phi),  cos(phi),   -cos(lambda)*sin(phi);
        cos(lambda),            0,          -sin(lambda)];

end