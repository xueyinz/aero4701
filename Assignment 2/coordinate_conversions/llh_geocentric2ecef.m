%% 440305585
% AERO4701
%
% Convert geocentric latitude, longitude, height coordinates to
% Earth-Centred Earth-Fixed coordinates
%
% inputs:	pos_LLH = [lat; long; alt] = LLH geocentric coordinates [rad, rad, m]
% outputs:  pos_ECEF = [x; y; z] = ECEF coordinates [m, m, m]

function pos_ECEF = llh_geocentric2ecef(pos_LLH)

    r_Earth = 6378137;
    
    lambda = pos_LLH(1);
    phi = pos_LLH(2);
    R = pos_LLH(3) + r_Earth;
    
    pos_ECEF = [R*cos(lambda)*cos(phi);
        R*cos(lambda)*sin(phi);
        R*sin(lambda)];

end