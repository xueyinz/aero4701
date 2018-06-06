%% 440305585
% AERO4701
% Assignment 3
%
% get_energy_ellipsoid.m

function ellipsoid_surf = get_energy_ellipsoid(E_total, I)

    % mesh size for the ellipsoid surf
    global surf_size;

    % semi-axis lengths
    ellipsoid_surf.a = sqrt(2 * E_total ./ I.xx);
    ellipsoid_surf.b = sqrt(2 * E_total ./ I.yy);
    ellipsoid_surf.c = sqrt(2 * E_total ./ I.zz);
    
    % ellipsoid surface
    [ellipsoid_surf.x, ellipsoid_surf.y, ellipsoid_surf.z] = ...
        ellipsoid(0, 0, 0, ellipsoid_surf.a, ellipsoid_surf.b, ellipsoid_surf.c, surf_size);

end