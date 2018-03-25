function [xpk1, Pk1, xk1] = UKF(esti, xhatk, Pk, y, xk, b, s, dt)

% INPUT:
% Esti - Estimation parameters
% y - Measurements
% qk - Last estimated quaternion
% b - Magnetic field model vector
% Pk - Last estimated covariance matrix
% dt - Time step
% xhatk - Last estimated USQUE state

% OUTPUT:
% xhatk1 - New estimated USQUE state
% Pk1 - New estimated covariance matrix
% qk1 - New estimated quaternion
% x - Estimated attitude state

%% Initialise

% Extract the components from the state matrix
omegaM = y(1:3);
if all(isnan(y(10:12)))
    M = y(4:9);
    noSun = true;
else
     M = y(4:12);
    noSun = false;
end
qk = xk(4:7);

% Scaling parameter f to go to generalised Rodrigues parameters
a = esti.a; lambda = esti.lambda;
f = (a+1)*2;

% Find the variances
var = esti.var;

% Calcualte Q
Q = [(var.v-var.u*dt^2/6)*eye(3) zeros(3); zeros(3) var.u*eye(3)]*dt/2;

% Calculate the sigma points
n = size(xhatk,1); NSig = 2*n+1;
sigk = chol((n+lambda)*(Pk + Q))';

% Create the columns of the chi matrix from the sigma points, as each sigma
% column plus the estimated state post update k
chiSigk = xhatk*ones(1,NSig)+[zeros(n,1), sigk, -sigk];

% Preallocate the next chi matrix and mean observation
chiSigk1 = zeros(6,NSig);
ySigk1 = zeros(numel(M),NSig);

% Loop through each sigma point
for idx = 1:NSig
    
    % Find the error quaternion from the generalised rodriques parameter 
    dq4k = (-a*norm(chiSigk(1:3,idx))^2 + ...
        f*sqrt(f^2+(1-a^2)*norm(chiSigk(1:3,idx))^2))/...
        (f^2 + norm(chiSigk(1:3,idx))^2);
    
    % Find the rest of the error quaternion
    dq13k = (a+dq4k)*chiSigk(1:3,idx)/f;
    
    % Find the quaternion error
    qki = quatmult([dq13k; dq4k], qk);
    
    % Find the estimated angular velocities
    omegaki = omegaM - chiSigk(4:6,idx);
    
    % Angular velocity in orbit to control frame
    woc = omegaki - Rob(qki)*[0; esti.wo; 0];
    
    % Evaluate the OMEGA matrix
    psik = (sin(0.5*norm(woc)*dt)/norm(woc))*woc;
    OM = [cos(0.5*dt*norm(woc))*eye(3) - crossPr(psik), psik;
        -psik', cos(0.5*dt*norm(woc))];

    % Kinematics
    qk1mi = OM*qki;
    
    % Propogate the error quaternions
    if idx == 1
        qkm1 = qk1mi; 
    end
    dqk1mi = quatmult(qk1mi,quatinv(qkm1));
    
    % Propogate the sigma points
    chiSigk1(1:3,idx) = f*dqk1mi(1:3)/(a+dqk1mi(4));
    chiSigk1(4:6,idx) = chiSigk(4:6,idx);
    
    % Calcualte the observation points
    if noSun
        ySigk1(:,idx) = [Rob(qk1mi)*b; Rob(qk1mi)*b];
    else
        ySigk1(:,idx) = [Rob(qk1mi)*b; Rob(qk1mi)*b; Rob(qk1mi)*s];
    end
end

% Calculate the predicted mean, sum using the W matrix multiplication
W = ones(NSig,1)/(2*(n+lambda));
W(1,1) = lambda/(n+lambda);
xmk1 = chiSigk1*W;
ymk1 = ySigk1*W;

% Find the innovation covariance
Rk1 = var.b*eye(6); 

% Find predicted covariance
Pmk1 = Q;
Pvv = Rk1;
Pxy = zeros(n,6);
for idx = 1:NSig
    xdiff = (chiSigk1(:,idx) - xmk1);
    ydiff = (ySigk1(:,idx) - ymk1);
    Pmk1 = Pmk1 + xdiff*xdiff'*W(idx,1);
    Pvv = Pvv + ydiff*ydiff'*W(idx,1);
    Pxy = Pxy + xdiff*ydiff'*W(idx,1);
end

% Find the gain matrix
K = Pxy/Pvv;

% Update
Pk1 = Pmk1 - K*Pxy';
xpk1 = xmk1 + K*(M-ymk1);

% Calcualte the new quaternion
dq4k1 = (-a*norm(xpk1(1:3))^2 + ...
    f*sqrt(f^2+(1-a^2)*norm(xpk1(1:3))^2))/...
    (f^2 + norm(xpk1(1:3))^2);

% Find the rest of the error quaternion
dq13k1 = (a+dq4k1)*xpk1(1:3)/f;

% Find the quaternions
qk1 = quatmult([dq13k1; dq4k1], qkm1);
qk1 = qk1/norm(qk1);

% Find the state output
xk1 = [omegaM-xpk1(4:6,1)-Rob(qk1)*[0; esti.wo; 0]; qk1];

% Reset xhat to be zeros again
xpk1(1:3) = zeros(3,1);

function n = quatmult(a,b)

% Calculate the components
n = [a(4)*b(1) + a(1)*b(4) + a(2)*b(3) - a(3)*b(2);
    a(4)*b(2) - a(1)*b(3) + a(2)*b(4) + a(3)*b(1);
    a(4)*b(3) + a(1)*b(2) - a(2)*b(1) + a(3)*b(4);
    a(4)*b(4) - a(1)*b(1) - a(2)*b(2) - a(3)*b(3)];

function q = quatinv(q)

% Find the inverse
q = [-q(1:3); q(4)];

function A = crossPr(a)

% Find the cross product
A = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];

function R = Rob(q)

% Rotation matrix
R = [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2, ...
    2*(q(1)*q(2) + q(3)*q(4)), 2*(q(1)*q(3) - q(2)*q(4));
    2*(q(1)*q(2) - q(3)*q(4)), -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,...
    2*(q(2)*q(3) + q(1)*q(4));
    2*(q(1)*q(3) + q(2)*q(4)), 2*(q(2)*q(3) - q(1)*q(4)), ...
    -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2];
