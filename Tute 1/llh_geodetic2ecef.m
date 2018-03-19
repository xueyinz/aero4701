%% 440305585
% AERO4701
%
% Convert geodetic latitude, longitude, height coordinates to
% Earth-Centred Earth-Fixed coordinates
%
% inputs:	pos_LLH = [lat; long; alt] = LLH geodetic coordinates [rad, rad, m]
% outputs:  pos_ECEF = [x; y; z] = ECEF coordinates [m, m, m]

function pos_ECEF = llh_geodetic2ecef(pos_LLH)

    if ((size(pos_LLH,1) ~= 3) || (size(pos_LLH,2) ~= 1))
        error('Check that dimensions of inputs to llh_geodetic2ecef match: ([lat;long;alt])')
    end

    global r_Earth;
    global e;
    
    lambda = pos_LLH(1);
    phi = pos_LLH(2);
    h = pos_LLH(3);
    N = r_Earth/sqrt(1 - (e^2)*(sin(lambda))^2);
    
    pos_ECEF = [(N+h)*cos(lambda)*cos(phi);
        (N+h)*cos(lambda)*sin(phi);
        N*(1 - e^2)*sin(lambda)];

end