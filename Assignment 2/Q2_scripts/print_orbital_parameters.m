%% 440305585
% AERO4701
%
% Print out the orbit parameters to the command window.
%
% inputs:   satellite = struct with parameters:
%                           a = semi-major axis
%                           e = eccentricity
%                           i = inclination
%                           Omega = right ascension of ascending node
%                           w = argument of perigee
%                           M0 = mean anomaly at epoch
%                           p = semilatus rectum
%                           n = mean motion
%           type = where the orbital parameters come from (string)
% outputs:  command window message

function print_orbital_parameters(original_params, initial_params, refined_params)

    struct_array = [original_params, initial_params, refined_params];
    table = struct2table(struct_array);
    disp(table);
    fprintf('^ Units for the orbital parameters (tabulated above) are in metres and radians.\n');

end