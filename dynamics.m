function [dS] = dynamics(~, State, J)
% This function represents both the rotational dynamics model, as well as
% the quaternion (attitude) kinematics model.  This is the entire system of
% physics to be simulated:

    % Recover State Variables:
    q = State(1:4);
    omega = State(5:7);
    
    % Quaternion (Attitude) Kinematics:
    Bq = zeros(4,3);
    Bq(1:3,:) = CrossProdMat(q(1:3)) + diag([q(4), q(4), q(4)]);
    Bq(4,:) = -q(1:3);
    dq = (1/2)*Bq*omega;

    % Rotational Dynamics (Euler's Equations):
    L_applied = zeros(3,1);
    L_dist = zeros(3,1);
    L = L_dist + L_applied;
    
    dOmega = J\(L - cross(omega,(J*omega)));

    dS = [dq; dOmega];
end
