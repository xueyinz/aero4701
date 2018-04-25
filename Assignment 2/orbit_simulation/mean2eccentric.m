%% 440305585
% AERO4701
%
% Numerically solve for eccentric anomaly E until the error between the
% current and new estimates for E converges below the accepted threshold.
%
% inputs:   M = mean anomaly [rad]
%           e = eccentricity
% outputs:  E = eccentric anomaly [rad]


function E_out = mean2eccentric(M_in, e)

    global E_threshold;

    E_out = zeros(1, length(M_in));
    for ii = 1:length(M_in)

        M = M_in(ii);
        E = M;

        temp = E - e*sin(E) - M;

        while abs(temp) > E_threshold
            E = E - (E - e*sin(E)- M)/(1-e*cos(E));
            temp = E - e*sin(E) - M;
        end
        
        E_out(ii) = E;
        
    end

end