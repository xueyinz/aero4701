%% 440305585
% AERO4701
%
% Generate pseudorandom error for a given deviation amount.
%
% input:    deviation = error deviation amount about 0
% output:   error = random error about 0

function noisy_measurements = add_random_error(measurements, deviation)

    error = -deviation + rand(size(measurements))*2*deviation;
    noisy_measurements = measurements + error;
    
end