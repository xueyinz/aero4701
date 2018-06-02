function out = dcm(in)
% Direction cosine matrix from quaternions
% The function takes in a [3x1] (psi,theta, phi) 
% Euler angle vector and outputs a [3x3] DCM
% matrix converting from Earth axis to body

psi = in(1); theta = in(2); phi = in(3); 

t1 = cos(psi)*cos(theta);
t2 = sin(psi)*cos(theta);
t3= -sin(theta);

m1 = cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi);
m2 = sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi);
m3 = cos(theta)*sin(phi);

b1 = cos(psi)*sin(theta)*cos(phi)+sin(psi)*sin(phi);
b2 = sin(psi)*sin(theta)*cos(phi)- cos(psi)*sin(phi);
b3 = cos(theta)*cos(phi);

out = [t1 t2 t3;m1 m2 m3;b1 b2 b3];