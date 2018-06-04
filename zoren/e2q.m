function Quats = E2Q(Eulers)
% Converts Euler angles to Quaternion components
% Quats = Euler2Quat(Eulers)

[~, n] = size(Eulers);

Quats = zeros(4, n);

for l = 1:n
    
    phi = Eulers(1, l);
    theta = Eulers(2, l);
    psi = Eulers(3, l);
    
    q0 = cos(psi/2)*cos(theta/2)*cos(phi/2) + sin(psi/2)*sin(theta/2)*sin(phi/2);
    q1 = cos(psi/2)*cos(theta/2)*sin(phi/2) - sin(psi/2)*sin(theta/2)*cos(phi/2);
    q2 = cos(psi/2)*sin(theta/2)*cos(phi/2) + sin(psi/2)*cos(theta/2)*sin(phi/2);
    q3 = -cos(psi/2)*sin(theta/2)*sin(phi/2) + sin(psi/2)*cos(theta/2)*cos(phi/2);
    
    Quats(:, l) = [q0; q1; q2; q3];
    
end