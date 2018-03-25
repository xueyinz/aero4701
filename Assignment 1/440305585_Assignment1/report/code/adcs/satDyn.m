function Xdot = satDyn(X,m,b,sat)

% Initialise the state rates
Xdot = zeros(13,1);

% Calculate the linear state rates
Xdot(1:3) = X(4:6);
[th,phi,r] = cart2sph(X(1),X(2),X(3));
g = -3.986004418e14/r^2;
[gx,gy,gz] = sph2cart(th,phi,g);
Xdot(4:6) = [gx; gy; gz];

% Cross product matrix
Omega = @(x) [0, x(3), -x(2), x(1);
    -x(3), 0, x(1), x(2);
    x(2), -x(1), 0, x(3);
    -x(1), -x(2), -x(3), 0];

% Find the rotation matrix
R = @(q) [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2, ...
    2*(q(1)*q(2) + q(3)*q(4)), 2*(q(1)*q(3) - q(2)*q(4));
    2*(q(1)*q(2) - q(3)*q(4)), -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,...
    2*(q(2)*q(3) + q(1)*q(4));
    2*(q(1)*q(3) + q(2)*q(4)), 2*(q(2)*q(3) - q(1)*q(4)), ...
    -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2];

% Angular velocity in orbit to control frame
woc = X(11:13) - R(X(7:10))*[0; sat.wo; 0];

% Kinematics
Xdot(7:10) = 0.5*Omega(woc)*X(7:10);

% Disturbance torque - gravity gradient
bko = R(X(7:10))*[0; 0; -1];
Tgg = 3*sat.wo^2*cross(bko,sat.I*bko); Tdist = Tgg;

% Atmospheric drag
drag = sat.drag;
if drag.on
    
    % Find the velocity vector in body frame
    vuv = R(X(7:10))*[1; 0; 0];
    v = norm(X(4:6));
    
    % Face normals
    normals = [1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1];
    
    % Shadowing function
    vdn = dot(normals,repmat(vuv,1,6))';
    vdn(vdn < 0) = 0;
    
    % Magnitude of drag
    f = -0.5*drag.rho*v^2*drag.Cd*repelem(drag.faces,2,1).*vdn;
    
    % Find the drag on each face
    Td = 0;
    for idx = 1:6
      Td = Td + cross(drag.rcp(:,idx),vuv*f(idx));
    end

    % Langmuir probes
    % Normal is the projection of vuv onto the plane normal to [0; 0; 1]
    normals = cross(drag.langVec,cross(vuv,drag.langVec)); 
    normals = normals/norm(normals);
    
    % Angle of the probes wrt the velocity vector
    vdn = dot(normals,vuv);
    
    % Magnitude of the force
    f = -0.5*drag.rho*v^2*drag.CdL*drag.langmuirS*vdn;
    
    % Find total torque
    Td = Td + cross(drag.rcpL(:,1),vuv*f)...
        + cross(drag.rcpL(:,2),vuv*f);
    Tdist = Tdist + Td;
end

% Torque from the coils
Tcoil = cross(m,b);

% Dynamics
Xdot(11:13) = sat.I\(-cross(X(11:13),sat.I*X(11:13)) + Tdist + Tcoil);

