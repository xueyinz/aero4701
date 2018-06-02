function QuatDot = quatRates(Quats, BodyRates)
% Evaluate Quaternion rates using initial Quaternion and body rates
% QuatDot = QuatRates(Quats, BodyRates)

[~, n] = size(Quats);

QuatDot = zeros(4, n);

for l = 1:n
    
    q0 = Quats(1, l);
    q1 = Quats(2, l);
    q2 = Quats(3, l);
    q3 = Quats(4, l);
        
    temp = 1/2*[-q1, -q2, -q3;...
        q0, -q3, q2;...
        q3, q0, -q1;...
        -q2, q1, q0];
    
    QuatDot(:, l) = temp*BodyRates(:, l);
    
end