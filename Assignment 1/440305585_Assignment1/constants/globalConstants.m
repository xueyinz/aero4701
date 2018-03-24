%% 440305585
% AERO4701
%
% Constants

global omega;
omega = 7.292115e-5;        % Earth rotation rate [rad/sec]

global r_Earth;
r_Earth = 6378137;          % Earth radius [m]

global a;
a = 6378137;                % Earth semimajor axis [m]

global e;
e = 0.08181919;             % Earth eccentricity

global mu_Earth;
mu_Earth = 3.986004418e14;  % Earth standard gravitational parameter [m3.s-2]

global sec_per_day;
sec_per_day = 24*60*60;     % seconds per day [s]

global E_threshold;
E_threshold = 0.0001;       % error threshold when solving for E

global dt;
dt = 10;                    % time steps [s]

global t_VE;
t_VE = [2018 3 20 16 16 0]; % time since last Vernal Equinox alignment (nearest min) [year month day hour minute seconds]

