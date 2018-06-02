function xdot = f(x, u)
    

    global I_s; 
    global I_w; 
    
    % Quaternions
    q = x(1:4);
    
    % Spacecraft angular velocity
    w_s = x(5:7); 
    
    % Wheel velocity
    w_w = x(8:10); 
    
    % Wheel acceleration 
    w_w_dot = u; 
    xdot = NaN(10, 1); 
    
    % Quaternion rates eq. 
    xdot(1:4) = 0.5*quatmultiply(q', [0; w_s]')'; 
    % Euler equation
    xdot(5:7) = I_s\(cross(I_s*w_s + I_w*w_w, w_s) - I_w*w_w_dot); 
    % Update RW ang. vel. 
    xdot(8:10) = w_w_dot; 
    
end