%% 440305585
% AERO4701
%
% Generate pseudorandom error for a given deviation amount.
%
% input:    deviation = error deviation amount about 0
% output:   error = random error about 0

function error = random_error(deviation)

    error = -deviation + rand()*2*deviation;

end