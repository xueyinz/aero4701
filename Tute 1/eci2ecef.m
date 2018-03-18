function ECEF = eci2ecef(pos_ECI, t)

% convert Earth-Centred Inertial coordinate to Earth-Centred Earth-Fixed
% coordinate

omega = 7.292115e-5;
C = [cos(omega*t) sin(omega*t) 0; -sin(omega*t) cos(omega*t) 0; 0 0 1];
ECEF = C*pos_ECI;