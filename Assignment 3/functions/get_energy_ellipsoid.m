%% 440305585
% AERO4701
% Assignment 3
%
% get_energy_ellipsoid.m

function energy_ellipsoid = get_energy_ellipsoid(E_total, I)

    % mesh size for the ellipsoid surf
    global surf_size;

    % semi-axis lengths
    energy_ellipsoid.a = sqrt(2 * E_total ./ I.xx);
    energy_ellipsoid.b = sqrt(2 * E_total ./ I.yy);
    energy_ellipsoid.c = sqrt(2 * E_total ./ I.zz);
    
    % ellipsoid surface
    [energy_ellipsoid.x, energy_ellipsoid.y, energy_ellipsoid.z] = ...
        ellipsoid(0, 0, 0, energy_ellipsoid.a, energy_ellipsoid.b, energy_ellipsoid.c, surf_size);

end