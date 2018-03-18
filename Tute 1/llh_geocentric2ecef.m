function ECEF = llh_geocentric2ecef(LLH)
r_Earth = 6378137;
lambda = deg2rad(LLH(1));
phi = deg2rad(LLH(2));
R = LLH(3) + r_Earth;
ECEF = [R*cos(lambda)*cos(phi); R*cos(lambda)*sin(phi); R*sin(lambda)];