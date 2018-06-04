function X = Integrate(X_i, I, dt)


X1 = X_i;

k1 = stateRates(X1, I);
A = X1 + k1.*dt./2;

k2 = stateRates(A, I);
B = X1 + k2.*dt./2;

k3 = stateRates(B, I);
C = X1 + k3.*dt;

k4 = stateRates(C, I);

Xtemp = X1 + (k1 + 2.*k2 + 2.*k3 + k4).*dt./6;

% mu = sqrt(sum(Xtemp(4:7).^2));
mu = 1;

X(1:3) = Xtemp(1:3);
X(4:7) = Xtemp(4:7)./mu;


