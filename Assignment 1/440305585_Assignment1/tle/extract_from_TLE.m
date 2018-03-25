%% 440305585
% AERO4701
%
% Reads a TLE code from a text file and extracts the relevant satellite
% parameters
%
% inputs:   filename with TLE code inside
% outputs:  object with satellite parameters

function sat = extract_from_TLE(filename)

    global sec_per_day;
    global mu_Earth;

    fid = fopen(filename);

    %% read the first three lines of the text file
    sat.name = fgetl(fid);
    tle_line1 = fgetl(fid);
    tle_line2 = fgetl(fid);
    sat.number = tle_line1(3:7);
    
    %% date and time
    Epoch_year = str2double(tle_line1(19:20));
    Epoch_day = str2double(tle_line1(21:32));
    sat = date_time_TLE(sat, Epoch_year, Epoch_day);
    sat.t0 = [sat.year, sat.month, sat.day, sat.hour, sat.min, sat.sec];
    
    %% orbital parameters
    sat.i = deg2rad(str2double(tle_line2(9:16)));               % inclination [rad]
    sat.Omega = deg2rad(str2double(tle_line2(18:25)));          % right ascension of ascending node [rad]
    sat.e = str2double(['.' tle_line2(27:33)]);                 % eccentricity
    sat.w = deg2rad(str2double(tle_line2(35:42)));              % ascending node/argument of perigee (omega) [rad]
    sat.M0 = deg2rad(str2double(tle_line2(44:51)));             % mean anomaly at Epoch [rad]
    n = str2double(tle_line2(53:63));                           % mean motion [rev/day]
    sat.n = n*2*pi/sec_per_day;                                 % mean motion [rad/s] 
    
    %% other useful parameters
    sat.a = (mu_Earth/(sat.n^2))^(1/3);                         % semi-major axis [m]
    sat.p = sat.a*(1 - (sat.e)^2);                              % semilatus rectum [m]
    sat.period = sec_per_day/n;                                 % orbital period [s]
    
end