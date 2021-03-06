%% 440305585
% AERO4701
%
% herrick_gibbs.m
%
% Use the Herrick-Gibbs technique to get the initial orbit determination,
% i.e. use three consecutive position measurements to determine the
% velocity at the second measurement.
%
% inputs:   pos_ECI = [x; y; z] = satellite position measurements in ECI frame [m; m; m]
%           t = time vector corresponding to pos_ECI
% outputs:  v_estimates = array estimates in ECI frame [m; m; m; m/s; m/s; m/s]
%                       = there are two less estimates than there are measurements
%                       = there are 3 rows for position, 3 rows for velocity

function initial_estimates = herrick_gibbs(pos_ECI, t)

    global mu_Earth
    
    n_estimates = size(pos_ECI,2) - 2;          % subtract 2 because the first and last measurements don't have velocity estimates
    initial_estimates = NaN(6, n_estimates);    % array of estimate results
    
    g = zeros(1,3);
    h = zeros(1,3);
    h_gain = mu_Earth/12;
    for ii = 1:n_estimates
        
        r = pos_ECI(:, ii:ii+2);            % 3x3 matrix
        if sum(isnan(t(ii:ii+2))) > 0
            continue
        end
        r_norm = rssq(r, 1);                % 1x3 matrix        
        
        t_12 = t(ii+1) - t(ii);
        t_13 = t(ii+2) - t(ii);
        t_23 = t(ii+2) - t(ii+1);
        
        g(1) = t_23/(t_12*t_13);
        g(3) = t_12/(t_23*t_13);
        g(2) = g(1) - g(3);
        
        h(1) = t_23;
        h(3) = t_12;
        h(2) = h(1) - h(3);
        
        d = g + h*h_gain./(r_norm.^3);
        d(1) = -d(1);

        initial_estimates(:,ii+1) = [r(:,2); sum(repmat(d,3,1).*r, 2)];
        
    end

end