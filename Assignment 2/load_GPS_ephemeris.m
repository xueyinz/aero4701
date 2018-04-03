%% 440305585
% AERO4701
%
% input: [Sat. #, a (m), e, i (deg), Omega (deg), w (deg), M0 (deg), t0 (s)]

function GPS_ephemeris = load_GPS_ephemeris(text_file)

    global mu_Earth;

    fid = fopen(text_file);
    GPS_ephemeris = [];
    
    while ~feof(fid)
        tline = fgetl(fid);
        tline = str2double(strsplit(tline, '  '));
        tline = tline(:, 2:end);            % get rid of empty column
        sat.number = tline(1);              % satellite number
        sat.a = tline(2);                   % satellite semi-major axis [m]
        sat.e = tline(3);                   % satellite eccentricity
        sat.i = deg2rad(tline(4));          % satellite inclination [rad]
        sat.Omega = deg2rad(tline(5));      % satellite right ascension of ascending node [rad]
        sat.w = deg2rad(tline(6));          % satellite argument of perigee [rad]
        sat.M0 = deg2rad(tline(7));         % satellite mean anomaly at Epoch [rad]
        sat.t0 = tline(8);                  % satellite Epoch time [s]
        sat.n = sqrt(mu_Earth/(sat.a^3));   % satellite mean motion [rad/s]
        sat.p = sat.a*(1 - (sat.e)^2);                              % semilatus rectum [m]
        GPS_ephemeris = [GPS_ephemeris; sat];
    end
    
end