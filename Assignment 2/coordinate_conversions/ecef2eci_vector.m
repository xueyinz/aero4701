%% 440305585
% AERO4701
%
% Convert Earth-Centred Earth-Fixed coordinates into Earth-Centred Inertial
% coordinates
%
% inputs:   pos_ECEF = [x; y; z] = ECEF coordinates [m, m, m]
%           t = time since vernal equinox
% outputs:  pos_ECI = [x; y; z] = ECI coordinates [m, m, m]

function pos_ECI = ecef2eci_vector(pos_ECEF, t)

    global omega;
    
    cw = cos(omega*t);
    sw = sin(omega*t);
    
    x = pos_ECEF(1,:).*cw - pos_ECEF(2,:).*sw;
    y = pos_ECEF(1,:).*sw + pos_ECEF(2,:).*cw;
    
    pos_ECI = [x; y; pos_ECEF(3,:)];
    
end