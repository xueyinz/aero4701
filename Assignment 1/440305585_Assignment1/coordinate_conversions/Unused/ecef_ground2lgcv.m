%% 440305585
% AERO4701
%
% Convert Earth-Centred Earth-Fixed coordinates of a satellite to local 
% geocentric vertical coordinates relative to a ground station in latitude,
% longitude, and height.
%
% inputs:   pos_satellite_ECEF = [x; y; z] = ECEF coordinates [m, m, m]
%           pos_ground_LLH_geocentric = [lat; long; alt] = [rad, rad, m]
% outputs:  pos_LGCV = [x; y; z] = LGCV coordinates [m, m, m]

function pos_satellite_relative_LGCV = ecef_ground2lgcv(pos_satellite_ECEF, pos_ground_LLH_geocentric)

    if ((size(pos_satellite_ECEF,1) ~= 3) || (size(pos_satellite_ECEF,2) ~= 1) || (size(pos_ground_LLH_geocentric,1) ~= 3) || (size(pos_ground_LLH_geocentric,2) ~= 1))
        error('Check that dimensions of inputs to ecef_ground2lgcv match: ([x;y;z], [lat;long;alt])')
    end
    
    % convert groundstation geocentric LLH to ECEF
    pos_ground_LGCV = llh_geocentric2ecef(pos_ground_LLH_geocentric);
    pos_relative = pos_satellite_ECEF - pos_ground_LGCV;
    
    lambda = pos_ground_LLH_geocentric(1);      % latitude
    phi = pos_ground_LLH_geocentric(2);         % longitude
    
    % rotation matrix
    C = [-sin(lambda)*cos(phi), -sin(phi),  -cos(lambda)*cos(phi);
        -sin(lambda)*sin(phi),  cos(phi),   -cos(lambda)*sin(phi);
        cos(lambda),            0,          -sin(lambda)];
    
    pos_satellite_relative_LGCV = C\pos_relative;

end