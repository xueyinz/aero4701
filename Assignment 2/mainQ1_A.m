%% 440305585
% AERO4701
% Assignment 2
%
% mainQ1_A.m

orbits.pos_ECI = zeros(3, size(t_12hr,2));
orbits = repmat(orbits, 1, n_satellites);

% for each satellite, get the trajectory over a 12hr period
for ii = 1:n_satellites
    
    % current satellite
    sat = GPS_ephemeris(ii);
    
    % mean anomaly: t = t1 - t0
    M = sat.M0*ones_t + sat.n.*(t_12hr - sat.t0);
    
    % eccentric anomaly
    E = mean2eccentric(M, sat.e);
    
    % true anomaly
    theta = 2*atan2(sqrt((1 + sat.e)/(1 - sat.e))*sin(E/2), cos(E/2));
    
    % radius
    r = sat.p./(1 + sat.e.*cos(theta));
    
    % get ECI position vectors
    orbits(ii).pos_ECI = orbit2ECI(r, theta, sat.i, sat.Omega, sat.w);
    
end