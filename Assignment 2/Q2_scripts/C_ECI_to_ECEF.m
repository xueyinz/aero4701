%Conversion matrix from week 1-2 slides
%input is t_VE time since vernal equinox alignment
function C_out = C_ECI_to_ECEF(t_VE)
    %time since vernal equinox alignment
    t = t_VE;
    %rate of earth rotation (rad/s)
    w_ie = 7.292115e-5;
    %ECEF output from ECI input
    C_out = [cos(w_ie*t), sin(w_ie*t), 0;...
             -sin(w_ie*t), cos(w_ie*t), 0;...
              0          ,0           , 1];
end
