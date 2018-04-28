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
    
    % time since 
    t_since_epoch = t_12hr - sat.t0;
    
    % get ECI position vectors
    orbits(ii).pos_ECI = orbit2ECI(sat, t_since_epoch);
    
end