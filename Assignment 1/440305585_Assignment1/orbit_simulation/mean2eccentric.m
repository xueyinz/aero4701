%% 440305585
% AERO4701
%
% Numerically solve for eccentric anomaly E until the error between the
% current and new estimates for E converges below the accepted threshold.
%
% inputs:   M = mean anomaly [rad]
%           e = eccentricity
% outputs:  E = eccentric anomaly [rad]


function E = mean2eccentric(M, e)

    global E_threshold;

    E_old = M;
    E_new = E_old - (E_old-e*sin(E_old)- M)/(1 - e*cos(E_old));
    
    while(abs(E_new - E_old) > E_threshold)
        E_old = E_new;
		E_new = E_old - (E_old - e*sin(E_old) - M)/(1 - e*cos(E_old));
    end
    
    E = E_new;

end