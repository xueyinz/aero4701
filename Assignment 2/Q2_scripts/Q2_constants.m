%% 440305585
% AERO4701
%
% Q2_constants.m

global omega
omega = 7.292115e-5;                % Earth rotation rate [rad/sec]

global r_Earth
r_Earth = 6378137;                  % Earth radius [m]

global a
a = 6378137;                        % Earth semimajor axis [m]

global e
e = 0.08181919;                     % Earth eccentricity

global mu_Earth
mu_Earth = 3.986005e14;             % Earth standard gravitational parameter [m3.s-2]

global E_threshold;
E_threshold = 0.00001;              % error threshold when solving for eccentric anomaly

sat.a = 26400e3;                    % satellite semimajor axis [m]
sat.e = 0.001;                      % satellite eccentricity
sat.i = deg2rad(55);                % satellite inclination angle [rad]
sat.Omega = deg2rad(120);           % satellite right ascension of ascending node [rad]
sat.w = 0;                          % satellite argument of perigee [rad]
sat.M0 = 0;                         % satellite mean anomaly at epoch [rad]
sat.n = sqrt(mu_Earth/(sat.a^3));   % satellite mean motion [rad/s]
sat.p = sat.a*(1 - (sat.e)^2);      % satellite semilatus rectum [m]

sec_per_day = 24*60*60;             % seconds per day [s]
t_24hr = 0:100:sec_per_day;         % time vector for 24hr period

% ground station in LLH [rad; rad; m]
ground_LLH = [deg2rad(38.8086); deg2rad(-104.5257); 1901];

% index for elevation in polar coordinates
ELEVATION = 3;

% minimum elevation angle to be able to see above the horizon
min_elevation = deg2rad(7.5);

range_error = 5;                    % range tracking error deviation [m]
angle_error = deg2rad(0.0001);      % angle tracking error deviation [rad]