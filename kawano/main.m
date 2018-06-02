%  Daniel Kawano, Rose-Hulman Institute of Technology
%  Last modified:  Feb 02, 2016

clear all
close all
clc

%  Specify symbolic parameters and symbolic variables that are functions of
%  time:

syms omega1(t) omega2(t) omega3(t) psi(t) theta(t) phi(t)
syms lambda1 lambda2 lambda3
assume((lambda1 > 0) & (lambda2 > 0) & (lambda3 > 0))
  
%  (1)  Balance of angular momentum with respect to the T-handle's mass
%       center:

%  Evaluate the absolute time derivative of the angular momentum about the
%  mass center:

H = diag([lambda1, lambda2, lambda3])*[omega1; omega2; omega3];
omegaRF = [omega1; omega2; omega3];

DH = diff(H) + cross(omegaRF, H);

%  Sum moments about the mass center:

sumM = [0, 0, 0]';

%  Construct the ODEs for rotational motion of the T-handle:

ODEs1 = DH == sumM;

%  (2)  Angular velocity kinematics:

%  Relate the space-fixed basis {E1,E2,E3} to the principal corotational
%  basis {e1,e2,e3} using a 3-1-3 set of Euler angles:

R1 = [cos(psi), sin(psi), 0;        
      -sin(psi), cos(psi), 0;
      0, 0, 1]; 
  
R2 = [1, 0, 0;                    	
      0, cos(theta), sin(theta);
      0, -sin(theta), cos(theta)];

R3 = [cos(phi), sin(phi), 0;        
      -sin(phi), cos(phi), 0;
      0, 0, 1];  

%  Express the angular velocity vector in terms of {e1,e2,e3}:  
  
omega = (R3*R2*R1)*[0; 0; diff(psi)] + (R3*R2)*[diff(theta); 0; 0] + ...
        R3*[0; 0; diff(phi)];

%  Construct the ODEs for the orientation of the T-handle:

ODEs2 = [omega1; omega2; omega3] == omega;

%  (3)  Manipulate the system of ODEs into a form suitable for numerical 
%       integration:

%  The ODEs are already in first-order form, so simply compile the state 
%  equations and arrange the state variables:
     
StateEqns = simplify([ODEs1; ODEs2]);

StateVars = [omega1; omega2; omega3; psi; theta; phi]

%  Express the state equations in mass-matrix form, M(t,Y)*Y'(t) = F(t,Y):

[Msym, Fsym] = massMatrixForm(StateEqns, StateVars);

Msym = simplify(Msym)
Fsym = simplify(Fsym)

%  Convert M(t,Y) and F(t,Y) to symbolic function handles with the input 
%  parameters specified:

M = odeFunction(Msym, StateVars, lambda1, lambda2, lambda3);  
F = odeFunction(Fsym, StateVars, lambda1, lambda2, lambda3); 

%  Save M(t,Y) and F(t,Y):

save t_handle_ODEs.mat Msym Fsym StateVars M F