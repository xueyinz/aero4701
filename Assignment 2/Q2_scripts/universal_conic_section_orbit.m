function [Position_ECI, Velocity_ECI] = universal_conic_section_orbit(t,init_pos,init_vel)

% UNIVERSAL_CONIC_SECTION_ORBIT - determines the position and velocity of a
% satellite given the orbital parameters and an input time using a
% universal conic section solution to Kepler's equations (from initial
% position and velocity)
%
% Inputs: t - time since epoch in seconds
%         init_pos - Satellite position at epoch
%         init_vel - Satellite velocity at epoch
%
% Outputs: Position_ECI - ECI frame position at time t
%          Velocity_ECI - ECI frame velocity at time t
%
% See Week 6 Notes, Slide 23 and 24 for Details
%

% Constants
mu = 3.986005e14; % Gm (Earth's gravitational constant (m^3/s^2))
we = 7.2921151467e-5; % Earth's rotation rate

% Analytically solve change in eccentric anomoly equation for position and
% velocity at time t 
a_semi = 1/(2/norm(init_pos) - norm(init_vel)^2/mu);
phi = deltaE_solve(a_semi, mu, init_pos, init_vel, t);

% Compute intermediate variables
f = 1 - a_semi*(1 - cos(phi))/norm(init_pos);
g = t - (a_semi^(3/2))*(phi - sin(phi))/sqrt(mu);
r = a_semi*(1 - (1 - (norm(init_pos)/a_semi))*cos(phi)) + init_pos'*init_vel*sqrt(a_semi/mu)*sin(phi);
fdot = -sqrt(mu*a_semi)*sin(phi)/(r*norm(init_pos));
gdot = 1 - a_semi*(1 - cos(phi))/r;

% Compute position and velocity at time t
Position_ECI = f*init_pos + g*init_vel;
Velocity_ECI = fdot*init_pos + gdot*init_vel;

% end of universal_conic_section_orbit

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

