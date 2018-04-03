%% 440305585
% AERO4701
%
% Convert geocentric latitude, longitude, height coordinates to
% Earth-Centred Earth-Fixed coordinates
%
% inputs:	pos_LLH = [lat; long; alt] = LLH geocentric coordinates [rad, rad, m]
% outputs:  pos_ECEF = [x; y; z] = ECEF coordinates [m, m, m]

function pos_ECEF = llh_geocentric2ecef(pos_LLH)

    if ((size(pos_LLH,1) ~= 3) || (size(pos_LLH,2) ~= 1))
        error('Check that dimensions of inputs to llh_geocentric2ecef match: ([lat;long;alt])')
    end

    global r_Earth;
    
    lambda = pos_LLH(1);
    phi = pos_LLH(2);
    R = pos_LLH(3) + r_Earth;
    
    pos_ECEF = [R*cos(lambda)*cos(phi);
        R*cos(lambda)*sin(phi);
        R*sin(lambda)];

end