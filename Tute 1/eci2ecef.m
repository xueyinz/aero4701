%% 440305585
% AERO4701
%
% Convert Earth-Centred Inertial coordinates into Earth-Centred Earth-Fixed
% coordinates
%
% inputs:   pos_ECI = [x; y; z] = ECI coordinates [m, m, m]
%           t = time since vernal equinox
% outputs:  pos_ECEF = [x; y; z] = ECEF coordinates [m, m, m]

function pos_ECEF = eci2ecef(pos_ECI, t)

    if ((size(pos_ECI,1) ~= 3) || (size(pos_ECI,2) ~= 1) || (length(t) ~= 1))
        error('Check that dimensions of inputs to eci2ecef match: ([x;y;z], t)')
    end

    global omega;
    
    % rotation matrix
    C = [cos(omega*t),	sin(omega*t),   0;
        -sin(omega*t),  cos(omega*t),   0;
        0,              0,              1];
    
    pos_ECEF = C*pos_ECI;
    
end