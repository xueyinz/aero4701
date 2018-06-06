%% 440305585
% AERO4701
% Assignment 3
%
% Q1_answers.m

%% (c) maximum angular velocity

w_max.x = max(w(1).x);
w_max.y = max(w(1).y);
w_max.z = max(w(1).z);

%% (e) copper wire temperature

% angular momentum is conserved for the system
w_final = L(1).total(1)/max([I.xx I.yy I.zz]);

% rotational kinetic energy in the final stable state
E_minimum = (1/2) * (max([I.xx I.yy I.zz]) * w_final^2);

% loss in rotational kinetic energy
E_loss = E(1).total(1) - E_minimum;

% wire properties
wire.volume = pi * (wire.d/2)^2 * wire.l;
wire.mass = wire.density * wire.volume;             % [g]

% change in thermal energy for the wire
wire.delta_temp = E_loss / (wire.mass * wire.c);
%     wire.final_temp = wire.temp + wire.delta_temp;

%% results

fprintf('1. (a)\n');
fprintf('\tI_xx = %f kg.m^2\n', I.xx);
fprintf('\tI_yy = %f kg.m^2\n', I.yy);
fprintf('\tI_zz = %f kg.m^2\n', I.zz);

fprintf('1. (b)\n');
fprintf('\tL_total = %f kg.m^2/s\n', L(1).total(1));

fprintf('1. (c)\n');
fprintf('\tmax(w_x) = %f rad/s\n', w_max.x);
fprintf('\tmax(w_y) = %f rad/s\n', w_max.y);
fprintf('\tmax(w_z) = %f rad/s\n', w_max.z);

fprintf('1. (d)\n');
fprintf('\tE_rotational = %f J\n', E(1).total(1));

fprintf('1. (e)\n');
fprintf('\tChange in temperature of wire = %f celsius\n', wire.delta_temp);
