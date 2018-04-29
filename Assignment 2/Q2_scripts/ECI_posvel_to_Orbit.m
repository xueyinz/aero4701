% Uses Herrick's Gibbs to convert pos and vel into orbital parameters
% input = [X,Y,Z,Xdot,Ydot,Zdot] in ECI
% output = [a,e,i,Omega,omega,M] in m and degrees
function output = ECI_posvel_to_Orbit(input)
%Constants
G = 6.67408e-11 ;
massearth = 5.972e24;
mu = G * massearth;
rad2deg = 180/pi;
deg2rad = pi/180;


r2 = input(1:3);
r2dot = input(4:6);

%Semimajor axis (a)
ivec = [1,0,0];
kvec = [0,0,1];
a = inv(2/norm(r2) - norm(r2dot)^2/mu);

% Inclination of orbit, (inc)
W = cross(r2,r2dot)/norm(cross(r2,r2dot));
inc = acos(dot(W,kvec));

% Eccentricy (ecc)
e = 1/mu * ( ((norm(r2dot))^2 - mu/norm(r2))*r2 - dot(r2,r2dot) * r2dot);
ecc = norm(e);
N = cross(kvec,W)/norm(cross(kvec,W));

% Right ascension of the ascending node (Omega)
cos_Omega = dot(ivec,N);
sin_Omega = dot(cross(ivec,N),kvec);
% Make sure that the value is in a period of 360degrees
Omega = wrapTo2Pi(atan2(sin_Omega,cos_Omega));

% Argument of perigee (omega)
cos_omega = dot(N,e)/norm(e);
sin_omega = cross(N,e)/norm(e) * W;
% Make sure that the value is in a period of 360degrees
omega = wrapTo2Pi(atan2(sin_omega,cos_omega));

% True anomaly (theta)
cos_theta = dot(e,r2)/(norm(e)*norm(r2));
sin_theta = dot((cross(e,r2)/(norm(e) * norm(r2))),W);
% Make sure that the value is in a period of 360degrees
theta = wrapTo2Pi(atan2(sin_theta,cos_theta));

% Mean Anomaly (M)
E = 2 * atan(( (1 - ecc)/(1 + ecc) )^(0.5) * tan(theta/2) );
M = E - ecc * sin(E);

output = [a,ecc,inc,Omega,omega,M];
end
