%% 440305585
% AERO4701
%
% Gets transformation matrix to convert ECI into ECEF coordinates.
%
% inputs:   t = time since last vernal equinox passing [s]
% outputs:  C = 3x3 transformation matrix

function C = C_eci2ecef(t)

    global omega;

    C = [cos(omega*t),  sin(omega*t),	0;
         -sin(omega*t), cos(omega*t),	0;
         0          ,   0           ,   1];

end