function Eulers = Q2E(Quats)
% Converts quaternion components to Euler angles.
% Takes in a [4xN] quaternion matrix. Outputs a [3xN] Euler angle matrix.
% Eulers = Q2E(Quats)

[~, n] = size(Quats);

Eulers = zeros(3, n);

for l = 1:n
    
    q0 = Quats(1, l);
    q1 = Quats(2, l);
    q2 = Quats(3, l);
    q3 = Quats(4, l);
    
    phi = atan2((q2*q3 + q0*q1) , (q0^2 + q3^2 - 1/2));
    theta = atan2((q0*q2 - q1*q3) , sqrt((q0^2 + q1^2 - 1/2)^2 + ...
        (q1*q2 + q0*q3)^2));
    psi = atan2((q1*q2 + q0*q3) , (q0^2 + q1^2 - 1/2));
    
    Eulers(:, l) = [phi; theta; psi];
    
    
end