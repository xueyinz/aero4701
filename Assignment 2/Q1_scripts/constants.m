%% 440305585
% AERO4701
%
% constants.m

global omega;
omega = 7.292115e-5;                % Earth rotation rate [rad/sec]

global r_Earth;
r_Earth = 6378137;                  % Earth radius [m]

global a;
a = 6378137;                        % Earth semimajor axis [m]

global e;
e = 0.08181919;                     % Earth eccentricity

global mu_Earth;
mu_Earth = 3.986005e14;             % Earth standard gravitational parameter [m3.s-2]

global sec_per_day;
sec_per_day = 24*60*60;             % seconds per day [s]

global E_threshold;
E_threshold = 0.000001;             % error threshold when solving for eccentric anomaly

t_VE = 7347737.336;                 % time since last vernal equinox passage [s]

%% Q1.A

dt = 120;                           % time steps [s]
t_12hr = 0:dt:(sec_per_day/2);      % time vector for 12 hours

%% Q1.B

error_threshold = 1;                % error threshold for non-linear least squares method
max_iterations = 10;                % maximum number of iterations before stopping NLLS
XYZ = 1:3;                          % indices in state vector for x, y, z position
CB = 4;                             % index in state vector for clock bias

ground_LLH = [deg2rad(-34.76); deg2rad(150.03); 680];   % ground station LLH [deg, deg, m]

AZIMUTH = 2;                        % index constant for azimuth
ELEVATION = 3;                      % index constant for elevation
ALTITUDE = 3;                       % index constant for altitude

grey = 0.6;                         % for plotting purposes