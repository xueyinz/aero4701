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
