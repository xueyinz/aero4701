function LLH = ecef2llh(pos_ECEF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECEF2LLH - function converts a position vector in the Earth Centered
% Earth Fixed frame to latitude, longitude and altitude

% Parameters
% tol = 1e-6;
a = 6378137; % Earth semimajor axis (m)
e = 0.08181919; % Earth eccentricity
w = 7.292115e-5; % Earth rotation rate (rad/s)

x = pos_ECEF(1);
y = pos_ECEF(2);
z = pos_ECEF(3);

% initial values
h = 0;
N = a;
p = sqrt(x^2 + y^2);
lprev = 0;
lambda = 100;

% iterate to solve for lat,long,alt
for i = 1:100
    sinlambda = z/(N*(1 - e^2) + h);
    lambda = atan((z + (e^2)*N*sinlambda)/p);
    N = a/sqrt(1 - (e^2)*(sin(lambda)^2));
    Nprev = N;
    h = p/cos(lambda) - N;
end

% final value
LLH = [lambda; atan2(y,x); h];

% end of ecef2llh

