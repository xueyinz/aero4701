function Xdot = stateRates(X, inertiaMtr)

p = X(1);
q = X(2);
r = X(3);

q0 = X(4);
q1 = X(5);
q2 = X(6);
q3 = X(7);

Ixx = inertiaMtr(1);
Iyy = inertiaMtr(2);
Izz = inertiaMtr(3);

pDot = (Iyy-Izz)*q*r/Ixx;
qDot = (Izz-Ixx)*r*p/Iyy;
rDot = (Ixx-Iyy)*p*q/Izz;

qDotVec = quatRates([q0; q1; q2; q3], [p; q; r]);

Xdot = [pDot; qDot; rDot; qDotVec];

% % % p = X(1);
% % % q = X(2);
% % % r = X(3);
% % % 
% % % q0 = X(4);
% % % q1 = X(5);
% % % q2 = X(6);
% % % q3 = X(7);
% % % 
% % % Ixx = inertiaMtr(1);
% % % Iyy = inertiaMtr(2);
% % % Izz = inertiaMtr(3);
% % % 
% % % pDot = (Iyy-Izz)*q*r/Ixx;
% % % qDot = (Izz-Ixx)*r*p/Iyy;
% % % rDot = (Ixx-Iyy)*p*q/Izz;
% % % 
% % % qDotVec = quatRates([q0; q1; q2; q3], [p; q; r]);
% % % 
% % % Xdot = [pDot; qDot; rDot; qDotVec];

