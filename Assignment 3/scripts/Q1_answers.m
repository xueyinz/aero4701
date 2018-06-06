%% 440305585
% AERO4701
% Assignment 3
%
% Q1_answers.m

%% (e) copper wire temperature

% angular momentum is conserved for the system
w.final = L.total_initial/max([I.xx I.yy I.zz]);

% rotational kinetic energy in the final stable state
E.minimum = (1/2) * (max([I.xx I.yy I.zz]) * w.final^2);
% E.minimum = (1/2) * L.total_initial^2 / min([I.xx, I.yy, I.zz]);

% loss in rotational kinetic energy
E.loss = E.total(1) - E.minimum;

% wire properties
wire.volume = pi * (wire.d/2)^2 * wire.l;
wire.mass = wire.density * wire.volume;             % [g]

% change in thermal energy for the wire
wire.delta_temp = E.loss / (wire.mass * wire.c);
%     wire.final_temp = wire.temp + wire.delta_temp;

%% results

fprintf('1. (a)\n');
fprintf('\tI_xx = %f kg.m^2\n', I.xx);
fprintf('\tI_yy = %f kg.m^2\n', I.yy);
fprintf('\tI_zz = %f kg.m^2\n', I.zz);

fprintf('1. (b)\n');
fprintf('\tL_total = %f kg.m^2/s\n', L.total_initial);

fprintf('1. (c)\n');
fprintf('\tmax(w_x) = %f rad/s\n', w.x_max);
fprintf('\tmax(w_y) = %f rad/s\n', w.y_max);
fprintf('\tmax(w_z) = %f rad/s\n', w.z_max);

fprintf('1. (d)\n');
fprintf('\tE_rotational = %f J\n', E.total(1));

fprintf('1. (e)\n');
fprintf('\tChange in temperature of wire = %f celsius\n', wire.delta_temp);
