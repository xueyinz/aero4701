%% 440305585
% AERO4701
% Assignment 3
%
% Q3_polhode.

%% pre-allocate memory for arrays and structs

% mesh size for the ellipsoid surf
global surf_size;

% energy ellipsoid parameters
energy_ellipsoid.a = NaN;
energy_ellipsoid.b = NaN;
energy_ellipsoid.c = NaN;
energy_ellipsoid.x = NaN(surf_size + 1, surf_size + 1);
energy_ellipsoid.y = NaN(surf_size + 1, surf_size + 1);
energy_ellipsoid.z = NaN(surf_size + 1, surf_size + 1);
energy_ellipsoid = repmat(energy_ellipsoid, 1, 3);

momentum_ellipsoid.a = NaN;
momentum_ellipsoid.b = NaN;
momentum_ellipsoid.c = NaN;
momentum_ellipsoid.x = NaN(surf_size + 1, surf_size + 1);
momentum_ellipsoid.y = NaN(surf_size + 1, surf_size + 1);
momentum_ellipsoid.z = NaN(surf_size + 1, surf_size + 1);
momentum_ellipsoid = repmat(momentum_ellipsoid, 1, 3);

ellipsoid_limits = [0, 0, 0];

%% loop three times for rotation about the three axes

for ii = 1:3
    
    % get the semi-axis lengths for each ellipsoid
    energy_ellipsoid(ii) = get_energy_ellipsoid(E(ii).total(1), I);
    momentum_ellipsoid(ii) = get_momentum_ellipsoid(L(ii).total(1), I);
    
    % get limits for plotting
    ellipsoid_limits(ii) = max([energy_ellipsoid(ii).a, energy_ellipsoid(ii).b, energy_ellipsoid(ii).c ...
        momentum_ellipsoid(ii).a, momentum_ellipsoid(ii).b, momentum_ellipsoid(ii).c]);
    
end
