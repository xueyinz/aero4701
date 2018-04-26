%% 440305585
% AERO4701
%
% Convert local geocentric vertical coordinates relative to a ground 
% station in latitude, longitude, and height to Earth-Centred Earth-Fixed 
% coordinates of a satellite
%
% inputs:   pos_LGCV = [x; y; z] = LGCV coordinates [m, m, m]
%           pos_ground_LLH_geocentric = [lat; long; alt] = [rad, rad, m]
% outputs:  pos_satellite_ECEF = [x; y; z] = ECEF coordinates [m, m, m]

function pos_satellite_ECEF = lgcv2ecef_ground(pos_satellite_LGCV, pos_ground_LLH_geocentric)

    if ((size(pos_satellite_LGCV,1) ~= 3) || (size(pos_satellite_LGCV,2) ~= 1) || (size(pos_ground_LLH_geocentric,1) ~= 3) || (size(pos_ground_LLH_geocentric,2) ~= 1))
        error('Check that dimensions of inputs to ecef_ground2lgcv match: ([x;y;z], [lat;long;alt])')
    end
    
    lambda = pos_ground_LLH_geocentric(1);      % latitude
    phi = pos_ground_LLH_geocentric(2);         % longitude
    
    % rotation matrix
    C = [-sin(lambda)*cos(phi), -sin(phi),  -cos(lambda)*cos(phi);
        -sin(lambda)*sin(phi),  cos(phi),   -cos(lambda)*sin(phi);
        cos(lambda),            0,          -sin(lambda)];
    
    % convert groundstation geocentric LLH to ECEF
    pos_ground_LGCV = llh_geocentric2ecef(pos_ground_LLH_geocentric);
    pos_satellite_relative_LGCV = C*pos_satellite_LGCV;
    pos_satellite_ECEF = pos_ground_LGCV + pos_satellite_relative_LGCV;

end
