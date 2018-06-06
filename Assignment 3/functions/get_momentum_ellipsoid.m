%% 440305585
% AERO4701
% Assignment 3
%
% get_momentum_ellipsoid.m

function momentum_ellipsoid = get_momentum_ellipsoid(L_total, I)

    % mesh size for the ellipsoid surf
    global surf_size;

    % semi-axis lengths
    momentum_ellipsoid.a = L_total ./ I.xx;
    momentum_ellipsoid.b = L_total ./ I.yy;
    momentum_ellipsoid.c = L_total ./ I.zz;
    
    % ellipsoid surface
    [momentum_ellipsoid.x, momentum_ellipsoid.y, momentum_ellipsoid.z] = ...
        ellipsoid(0, 0, 0, momentum_ellipsoid.a, momentum_ellipsoid.b, momentum_ellipsoid.c, surf_size);

end