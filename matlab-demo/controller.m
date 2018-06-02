function u = controller(x)

    global I_s; 
    global I_w; 
    global K; 
    
    w_s = x(5:7); 
    
    % This controller can be shown to globally stabilize w_s. This is
    % because, if you consider V(w_s) = w_s' * I_s * w_s, the kinetic
    % energy, the time derivative of kinetic energy is 
    % dV/dt = - w_s * I_w *K * w_s, which is always negative as long as K
    % is a positive definite matrix (i.e. all eigenvalues are positive). 

    u = I_w \ cross(w_s, I_s*w_s) + K * w_s;  

    % Such a function, V, is called the Lyapunov function of the system,
    % which shows the stability. The condition is that, if V is a positive
    % definite function (V(0) = 0, but V(x) > 0 for all x != 0), and dV/dt
    % is negative definite (i.e. dV/dt < 0 for all x != 0) then x -> 0 as t -> inf.
    
end