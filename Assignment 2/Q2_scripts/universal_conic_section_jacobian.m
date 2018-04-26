function dxdxo = universal_conic_section_jacobian(pos, vel, init_pos, init_vel, delta_t)

% UNIVERSAL_CONIC_SECTION_JACOBIAN - Analytical solution to the derivative
% of the state transition function between the initial position and
% velocity and position and velocity and time t of a satellite in orbit
% (around the Earth) using the universal conic section method
%
% Inputs: init_pos - Initial ECI position at time to
%         init_vel - Initial ECI velocity at time to
%         pos - Final ECI position at time t = to + delta_t
%         vel - Final ECI velocity at time t = to + delta_t
%         pos - difference in time between t and to
%
% Outputs: dxdxo - Jacobian matrix of state transition function w.r.t
% initial position and velocity (6x6 matrix)
%
% See Week 6 Notes, Slide 29 for Details
%

% Constants
mu = 3.986005e14; % Gm (Earth's gravitational constant (m^3/s^2))

% Analytically solve change in eccentric anomoly equation for position and
% velocity at time t 
a_semi = 1/(2/norm(init_pos) - norm(init_vel)^2/mu);
phi = deltaE_solve(a_semi, mu, init_pos, init_vel, delta_t);

% Compute intermediate variables
f = 1 - a_semi*(1 - cos(phi))/norm(init_pos);
g = delta_t - (a_semi^(3/2))*(phi - sin(phi))/sqrt(mu);
r = a_semi*(1 - (1 - (norm(init_pos)/a_semi))*cos(phi)) + init_pos'*init_vel*sqrt(a_semi/mu)*sin(phi);
fdot = -sqrt(mu*a_semi)*sin(phi)/(r*norm(init_pos));
gdot = 1 - a_semi*(1 - cos(phi))/r;

% more intermediate variables
alpha = 1/a_semi;
chi = alpha*sqrt(mu)*delta_t + (pos'*vel - init_pos'*init_vel)/sqrt(mu);
u2 = (1 - cos(sqrt(alpha)*chi))/alpha;
u3 = (sqrt(alpha)*chi - sin(sqrt(alpha)*chi))/(alpha*sqrt(alpha));
u4 = (chi^2)/(2*alpha) - u2/alpha;
u5 = (chi^3)/(6*alpha) - u3/alpha;
c = (3*u5 - chi*u4 - sqrt(mu)*delta_t*u2)/sqrt(mu);

% State transition matrix elements
phi11 = (norm(pos)/mu)*(vel - init_vel)*(vel - init_vel)' + ...
    (norm(init_pos)^(-3))*(norm(init_pos)*(1 - f)*pos*init_pos' + c*vel*init_pos') + f*eye(3);
phi12 = (norm(init_pos)/mu)*(1 - f)*((pos - init_pos)*init_vel' - ...
    (vel - init_vel)*init_pos') + (c/mu)*vel*init_vel' + g*eye(3);
phi21 = -(norm(init_pos)^(-2))*(vel - init_vel)*init_pos' - (norm(pos)^(-2))*pos*(vel - init_vel)' - ...
    (mu*c/((norm(pos)^3)*(norm(init_pos)^3)))*pos*init_pos' ...
    + fdot*(eye(3) - (norm(init_pos)^(-2))*pos*pos' + (1/(mu*norm(pos)))*(pos*vel' - vel*pos')*pos*(vel - init_vel)');
phi22 = (norm(init_pos)/mu)*(vel - init_vel)*(vel - init_vel)' + ...
    (norm(init_pos)^(-3))*(norm(init_pos)*(1 - f)*pos*init_pos' - c*pos*vel') + gdot*eye(3);

% Compose STM elements
dxdxo = [phi11,phi12;phi21,phi22];

% end of universal_conic_section_jacobian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function phik = deltaE_solve(a,mu,r,rdot,delta_t)

% DELTAE_SOLVE - solve equation for change in eccentric anomoly using a
% Newton/Raphson method

tol = 1e-12;
phik = 0;
del_phi = 1;
i = 0;
while del_phi > tol
    phi_old = phik;
    B = (a^(3/2))/sqrt(mu);
    C = 1 - norm(r)/a;
    D = r'*rdot/sqrt(mu*a);
    fn = B*(phik - C*sin(phik) + D*(1 - cos(phik))) - delta_t;
    fndot = B - B*C*cos(phik) + B*D*sin(phik);
    phik = phik - (fn/fndot);
    del_phi = rem(abs(phi_old - phik), 2*pi);
    i = i + 1;
end

% end of deltaE_solve

