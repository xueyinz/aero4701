%% 440305585
% AERO4701
%
% Calculate the ECI position coordinates of the satellite using the orbital
% parameters inputted.
%
% inputs:   satellite = struct with parameters:
%                           e = eccentricity
%                           i = inclination
%                           Omega = right ascension of ascending node
%                           w = argument of perigee
%                           M0 = mean anomaly at epoch
%                           p = semilatus rectum
%                           n = mean motion
%           t_since_epoch = time elapsed since epoch t0 (coincides with M0)
% outputs:  pos_ECI = satellite position coordinates in ECI frame
%                       (3xN matrix where N is the number of time-steps)

function pos_ECI = orbit2ECI(satellit, t_since_epoch)

    % for code readability
    e = satellit.e;
    i = satellit.i;
    Omega = satellit.Omega;
    w = satellit.w;
    M0 = satellit.M0;
    p = satellit.p;
    n = satellit.n;
    
    % mean anomaly
    M = M0 + n.*(t_since_epoch);
    
    % eccentric anomaly
    E = mean2eccentric(M, e);

    % true anomaly
    theta = 2*atan2(sqrt((1 + e)/(1 - e))*sin(E/2), cos(E/2));

    % radius
    r = p./(1 + e.*cos(theta));
    
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