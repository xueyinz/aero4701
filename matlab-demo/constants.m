global I_s; 
global I_w; 
global K; 

% Spacecraft MoI
I_s = [2, 0.01, -0.01;
          0.01, 1, 0.02;
          -0.01, 0.02 6] * 10^-3; 

% Reaction wheel MoI
I_w = [5, 0, 0; 
            0, 5, 0;
            0, 0, 5] * 10^-6; 
        
 % Controller gain. 
K = 5 * eye(3) * 10^3; 
        
      
      