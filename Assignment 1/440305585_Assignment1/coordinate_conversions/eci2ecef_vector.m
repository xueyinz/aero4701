%% 440305585
% AERO4701
%
% Convert Earth-Centred Inertial coordinates into Earth-Centred Earth-Fixed
% coordinates
%
% inputs:   pos_ECI = [x; y; z] = ECI coordinates [m, m, m] (x, y, z are vectors)
%           t = time since vernal equinox alignment
% outputs:  pos_ECEF = [x; y; z] = ECEF coordinates [m, m, m]

function pos_ECEF = eci2ecef_vector(pos_ECI, t)

    global omega;

    x = pos_ECI(1,:).*cos(omega*t) + pos_ECI(2,:).*sin(omega*t);
    y = -pos_ECI(1,:).*sin(omega*t) + pos_ECI(2,:).*cos(omega*t);
    z = pos_ECI(3,:);
    
    pos_ECEF = [x; y; z];
    
end