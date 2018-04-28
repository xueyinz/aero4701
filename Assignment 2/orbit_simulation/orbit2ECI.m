%% 440305585
% AERO4701
%
% Calculate the ECI position coordinates of the satellite using the orbital
% parameters inputted.
%
% inputs:   r = radius vector [m]
%           theta = true anomaly vector [rad]
%           i = inclination angle scalar [rad]
%           Omega = right ascension of ascending node scalar [rad]
%           w = argument of perigee scalar [rad]
% outputs:  pos_ECI = satellite position coordinates in ECI frame
%                       (3xN matrix where N is the number of time-steps)

function pos_ECI = orbit2ECI(r, theta, i, Omega, w)
    
    % convert orbit into an Equatorial plane with Cartesian coordinates
    x = r.*cos(theta);
    y = r.*sin(theta);
    z = zeros(1, length(x));

    % transform the orbit according to its inclination, ascending node, and
    % argument of perigee
    C = [(cos(Omega)*cos(w)-sin(Omega)*sin(w)*cos(i)) (-cos(Omega)*sin(w)-sin(Omega)*cos(w)*cos(i)) (sin(Omega)*sin(i));
        (sin(Omega)*cos(w)+cos(Omega)*sin(w)*cos(i)) (-sin(Omega)*sin(w)+cos(Omega)*cos(w)*cos(i)) (-cos(Omega)*sin(i));
        (sin(w)*sin(i)) cos(w)*sin(i) cos(i)];      % transformation matrix
    pos_ECI = C*[x; y; z];

end