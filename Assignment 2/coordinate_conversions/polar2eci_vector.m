%% 440305585
% AERO4701
%
% Convert polar coordinates into ECI coordinates
%
% inputs:   pos_POLAR = [r, az, el] = polar coordinates [m, rad, rad]
%           pos_ground_LLH = ground station in [lat; long; alt] = [rad, rad, m]
%           t = time vector corresponding to pos_POLAR matrix
% outputs:	pos_ECI = [x; y; z] = ECI coordinates [m, m, m] (x, y, z are vectors)

function pos_ECI = polar2eci_vector(pos_POLAR, pos_ground_LLH, t)

    pos_LGCV = polar2cartesian_vector(pos_POLAR);
    pos_ECEF = lgcv_ground2ecef_vector(pos_LGCV, pos_ground_LLH);
    pos_ECI = ecef2eci_vector(pos_ECEF, t);
    
end