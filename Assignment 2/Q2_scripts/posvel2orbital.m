%% 440305585
% AERO4701
%
% Convert position and velocity estimates into orbital parameters.
%
% inputs:   posvel = satellite position and velocity [m; m; m; m/s; m/s; m/s]
% outputs:  orbital_parameters = [a e i Omega w M0]

function orbital_parameters = posvel2orbital(posvel)

    global mu_Earth
    
    r = posvel(1:3);                % position
    rdot = posvel(4:6);             % velocity
    rr = norm(r);
    V = norm(rdot);
    cross_rv = cross(r,rdot);
    unit_i = [1,0,0];               % unit vector in i direction
    unit_k = [0,0,1];               % unit vector in k direction
    
    % semi-major axis
    a = 1/(2/rr - V^2/mu_Earth);

    % eccentricity
    e_vector = 1/mu_Earth * ( (V^2 - mu_Earth/rr)*r - dot(r,rdot) * rdot);
    ee = norm(e_vector);
    
    % inclination
    W_vector = cross_rv/norm(cross_rv);
    i = acos(dot(W_vector,unit_k));

    % right ascension of the ascending node
    cross_kW = cross(unit_k,W_vector);
    N = cross_kW/norm(cross_kW);
    cos_Omega = dot(unit_i,N);
    sin_Omega = dot(cross(unit_i,N), unit_k);
    Omega = atan2(sin_Omega,cos_Omega);
    
    % argument of perigee
    cos_w = dot(N,e_vector)/ee;
    sin_w = cross(N,e_vector)/ee * W_vector;
    w = atan2(sin_w,cos_w);
    
    % true anomaly
    cos_theta = dot(e_vector,r)/(ee*rr);
    sin_theta = dot((cross(e_vector,r)/(ee * rr)), W_vector);
    theta = atan2(sin_theta,cos_theta);
    
    % mean anomaly
    E = 2*atan(sqrt((1 - ee)/(1 + ee))*tan(theta/2));
    M0 = E - ee * sin(E);

    orbital_parameters = [a ee i Omega w M0];
    
end