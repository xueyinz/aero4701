%% Space Engineering 3 Assignment 3
% Author: Zoren Liu
% SID: 440353416
%

clear; close all; clc;

%% Question 1

% Initial conditions
omega_x = 0;%pi/2*0.025; % Angular velocity (rad/s)
omega_y = 0;%pi/2*0.025;
omega_z = pi/2;


% Properties of rectangular plate
a = 0.230;   % Length (m)
b = 0.01;    % Height (m)
c = 0.11;   % Width (m)
m = 1;      % Mass (kg)

% a) Calculate moments of inertia
Ixx = m*(a^2+b^2)/12;   % Intermediate moment of inertia
Iyy = m*(a^2+c^2)/12;   % Maximum moment of inertia
Izz = m*(b^2+c^2)/12;   % Lowest moment of inertia

% b) Calculate magnitude of inertial angular momentum vector L_tot
Lx = Ixx*omega_x;
Ly = Iyy*omega_y;
Lz = Izz*omega_z;
L_tot = Lx + Ly + Lz;

% c) Calculate the maximum angular velocity in the three principle axes
% Since angular momentum is fixed, then the theoretical maximum angular
% velocity of an axis is when the angular velocities of the other two axes
% are zero
omega_x_temp = 0;
omega_y_temp = 0;
omega_z_temp = 0;

omega_x_max = (L_tot - Iyy*omega_y_temp - Izz*omega_z_temp)/(Ixx);
omega_y_max = (L_tot - Ixx*omega_x_temp - Izz*omega_z_temp)/(Iyy);
omega_z_max = (L_tot - Ixx*omega_x_temp - Iyy*omega_y_temp)/(Izz);

% d) Calculate the magnitude of the rotational kinetic energy for the three
% principal axes components in the body reference frame (Exx, Eyy, Ezz).
Exx = 1/2*Ixx*omega_x^2;
Eyy = 1/2*Iyy*omega_y^2;
Ezz = 1/2*Izz*omega_z^2;

% E = Exx + Eyy + Ezz;

% e)
% Copper wire
l_Cu = 17.7; % cm
d_Cu = 0.1; % cm
rho_Cu = 9; % g/cm^3
m_Cu = pi*(d_Cu/2)^2*l_Cu*rho_Cu;
temp_0 = 0; % Celcius
C_Cu = 0.385; % Joules/g


E_max = 1/2*L_tot*omega_z;
E_min = 1/2*L_tot*omega_y_max;
delta_E = E_max - E_min;

% % Assuming pure rotation in z-axis
% L_ini = Izz*omega_z;
% E_max = 1/2*L_ini*omega_z;
% omega_y_max_2 = L_ini/Iyy;
% E_min = 1/2*L_ini*omega_y_max_2;
% delta_E = E_max - E_min;

temp_end = delta_E/(C_Cu*m_Cu)+temp_0;

%% Question 2

plate.A = [-a/2 -c/2 -b/2];
plate.B = [a/2 -c/2 -b/2];
plate.C = [-a/2 c/2 -b/2];
plate.D = [-a/2 -c/2 b/2];
plate.E = [-a/2 c/2 b/2];
plate.F = [a/2 -c/2 b/2];
plate.G = [a/2 c/2 -b/2];
plate.H = [a/2 c/2 b/2];

origin = [0, 0, 0];

sim_type = 1;
axis_type = 1;
animStep = 100;
time = 60;
dt = 0.001;
tvec = (dt:dt:time)';

switch sim_type
    case 1
        P = [plate.A;plate.B;plate.F;plate.H;...
            plate.G;plate.C;plate.A;plate.D;...
            plate.E;plate.H;plate.F;plate.D;...
            plate.E;plate.C;plate.G;plate.B];
        plot3(P(:,1),P(:,2),P(:,3),'k', 'LineWidth', 2), hold on % original cube
        plot3(P(1:2,1),P(1:2,2),P(1:2,3),'b', 'LineWidth', 2)
        plot3(P(8:9,1),P(8:9,2),P(8:9,3),'r', 'LineWidth', 2)
        plot3(P(2:3,1),P(2:3,2),P(2:3,3),'g', 'LineWidth', 2); hold off
        %     case 2
        %         vert_p = [A;B;G;C;D;F;H;E];
        %         face_p = [1 2 3 4; 5 6 7 8; 1 2 6 5; 3 4 8 7; 1 4 8 5; 2 3 7 6];
        %         colour_p = [0 1 0; 0 1 0; 1 0 0; 1 0 0; 0 0 1; 0 0 1];
        %         P = patch('Vertices', vert_p, 'Faces', face_p,...
        %             'FaceVertexCData', colour_p, 'FaceColor', 'flat');
end
axis equal; grid on; axis([-0.20 0.20 -0.20 0.20 -0.20 0.20]); rotate3d on;

switch axis_type
    case 1 % About x
        omega_x = -pi/2;
        omega_y = 0.01;
        omega_z = 0.01;
    case 2 % About y
        omega_x = -0.1;
        omega_y = pi/2;
        omega_z = 0.1;
    case 3 % About z
        omega_x = -0.1;
        omega_y = 0.1;
        omega_z = pi/2;
end

ang_x = 0;
ang_y = 0;
ang_z = 0;

omega = zeros(length(tvec),3);
omega(1, :) = [omega_x, omega_y, omega_z];
ang(1, :) = [0, 0, 0];
counter = 0;

for n = 2:length(tvec)
    
    omega_x_dot = (Iyy-Izz)*omega(n-1, 2)*omega(n-1, 3)/Ixx;
    omega_y_dot = (Izz-Ixx)*omega(n-1, 3)*omega(n-1, 1)/Iyy;
    omega_z_dot = (Ixx-Iyy)*omega(n-1, 1)*omega(n-1, 2)/Izz;
    
    omega(n, 1) = omega(n-1, 1) + omega_x_dot*dt;
    omega(n, 2) = omega(n-1, 2) + omega_y_dot*dt;
    omega(n, 3) = omega(n-1, 3) + omega_z_dot*dt;
    
    ang(n, 1) = ang(n-1, 1) + omega(n, 1)*dt;
    ang(n, 2) = ang(n-1, 2) + omega(n, 2)*dt;
    ang(n, 3) = ang(n-1, 3) + omega(n, 3)*dt;
    
end

ang2(:, 1) = cumtrapz(omega(:, 1))*dt;
ang2(:, 2) = cumtrapz(omega(:, 2))*dt;
ang2(:, 3) = cumtrapz(omega(:, 3))*dt;

% Euler angles
dcmMtr = dcm([ang2(:, 2), ang2(:, 1), ang2(:, 3)]);
P_plot = zeros(size(P, 1), size(P, 2), length(dcmMtr));
for n = 1:length(tvec)
    P_plot(:, :, n) = P*dcmMtr(:, :, n);
end

% Quaternions
q = e2q([ang2(:, 2), ang2(:, 1), ang2(:, 3)]');
dcmMtr = dcmQuat(q');
P_plot = zeros(size(P, 1), size(P, 2), length(dcmMtr));
for n = 1:length(tvec)
    P_plot(:, :, n) = P*dcmMtr(:, :, n);
end

E(:, 1) = 1/2*Ixx*omega(:, 1).^2;
E(:, 2) = 1/2*Iyy*omega(:, 2).^2;
E(:, 3) = 1/2*Izz*omega(:, 3).^2;
Etot = E(:, 1) + E(:, 2) + E(:, 3);

ellips1(:, 1) = 2*Etot./Ixx;
ellips1(:, 2) = 2*Etot./Iyy;
ellips1(:, 3) = 2*Etot./Izz;

L(:, 1) = omega(:, 1).*Ixx;
L(:, 2) = omega(:, 2).*Iyy;
L(:, 3) = omega(:, 3).*Izz;
Ltot = sqrt(sum(L.^2, 2));

ellips2(:, 1) = Ltot.^2/Ixx^2;
ellips2(:, 2) = Ltot.^2/Iyy^2;
ellips2(:, 3) = Ltot.^2/Izz^2;

d = 2*Etot./Ltot.^2;

for o = 1:animStep:length(tvec)
%     o
    figure(1)
    switch sim_type
        case 1
            subplot(1,2,1)
            plot3(P_plot(:,1,o),P_plot(:,2,o),P_plot(:,3,o),'k', 'LineWidth', 2), hold on
            plot3(P_plot(1:2,1,o),P_plot(1:2,2,o),P_plot(1:2,3,o),'b', 'LineWidth', 2)
            plot3(P_plot(8:9,1,o),P_plot(8:9,2,o),P_plot(8:9,3,o),'r', 'LineWidth', 2)
            plot3(P_plot(2:3,1,o),P_plot(2:3,2,o),P_plot(2:3,3,o),'g', 'LineWidth', 2)
            plot3(0,0,0,'or')
            xlabel('y'); ylabel('x'); zlabel('z')
            axis equal; grid on; axis([-0.20 0.20 -0.20 0.20 -0.20 0.20]); rotate3d on; hold off
            drawnow
        case 2
            subplot(1,2,1)
            rotate(P, [0 1 0], rad2deg(omega(o, 1)*dt*animStep), origin);
            rotate(P, [1 0 0], rad2deg(omega(o, 2)*dt*animStep), origin);
            rotate(P, [0 0 1], rad2deg(omega(o, 3)*dt*animStep), origin);
            axis equal; grid on; axis([-0.20 0.20 -0.20 0.20 -0.20 0.20]); rotate3d on; hold off
            drawnow
    end
    
%     subplot(1,2,2)
%     plot(tvec(1:o), omega(1:o, 2), tvec(1:o), omega(1:o, 1), tvec(1:o), omega(1:o, 3),...
%         tvec(1:o), Ltot(1:o),'LineWidth', 1.5)
%     xlabel('Time (s)'); ylabel('\omega (rad/s)');
%     legend('\omega_x', '\omega_y', '\omega_z')
%     grid on
%     hold off
    
    
    
            [x1, y1, z1] = ellipsoid(0,0,0,sqrt(ellips1(o, 1)),sqrt(ellips1(o, 2)),sqrt(ellips1(o, 3)),30);
            [x2, y2, z2] = ellipsoid(0,0,0,sqrt(ellips2(o, 1)),sqrt(ellips2(o, 2)),sqrt(ellips2(o, 3)),30);
    
    %         [x1, y1, z1] = ellipsoid(0,0,0,Lsum(o),Lsum(o),Lsum(o),30);
    %         [x2, y2, z2] = ellipsoid(0,0,0,ellips2(o, 1),ellips2(o, 2),ellips2(o, 3),30);
    
            figure(2)
            surf(x1, y1, z1, 'edgecolor', 'none','FaceAlpha',0.3,'FaceColor',[0.9100 0.4100 0.1700])
             hold on;
            surf(x2, y2, z2, 'edgecolor', 'none','FaceAlpha',0.3,'FaceColor',[0 0 1])
            xlabel('\omega_y'); ylabel('\omega_x'); zlabel('\omega_z')
            axis([-5 5 -5 5 -10 10]); grid on;
    %         axis tight; axis equal; grid on;
            plot3(omega(:, 1),omega(:, 2),omega(:, 3),'r.')
    
            drawnow; hold off;
    
end



figure(2)
plot(tvec, omega(:, 1), tvec, omega(:, 2), tvec, omega(:, 3),...
    'LineWidth', 1.5)
xlabel('Time (s)'); ylabel('\omega (rad/s)');
legend('\omega_x', '\omega_y', '\omega_z')
grid on
figure(3)
plot(tvec, ang(:, 1), tvec, ang(:, 2), tvec, ang(:, 3),...
    'LineWidth', 1.5)
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('ang_x', 'ang_y', 'ang_z')
grid on



% figure(2)
% plot(tvec, rad2deg(omegaDB(:, 1)), tvec, rad2deg(omegaDB(:, 2)), tvec, rad2deg(omegaDB(:, 3)),...
%     'LineWidth', 1.5)
% xlabel('Time (s)'); ylabel('\omega (deg/s)');
% legend('\omega_x', '\omega_y', '\omega_z')
% grid on
% figure(3)
% plot(tvec, rad2deg(angDB(:, 1)), tvec, rad2deg(angDB(:, 2)), tvec, rad2deg(angDB(:, 3)),...
%     'LineWidth', 1.5)
% xlabel('Time (s)'); ylabel('Angle (deg)');
% legend('ang_x', 'ang_y', 'ang_z')
% grid on



