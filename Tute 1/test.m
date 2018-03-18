%% test script

% ECEF = eci2ecef([7000000; 6100000; 5100000], 180);
% disp(ECEF)

% LLH_geocentric = ecef2llh_geocentric([-4678515; 2593343; -3473810]);
% disp(LLH_geocentric)

% ECEF = llh_geocentric2ecef([-33; 151; 52]);
% disp(ECEF)

% LLH_geodetic = [-33; 151; 52];
% ECEF = llh_geodetic2ecef(LLH_geodetic);
% disp(ECEF)
% LLH_geodetic_new = ecef2llh_geodetic(ECEF);
% LLH_geodetic_new(1) = rad2deg(LLH_geodetic_new(1));
% LLH_geodetic_new(2) = rad2deg(LLH_geodetic_new(2));
% disp(LLH_geodetic_new)

% CARTESIAN = polar2cartesian([1000; -140; 60]);
% disp(CARTESIAN)

% POLAR = cartesian2polar([-383.02; -321.39; -866.03]);
% disp(POLAR)

