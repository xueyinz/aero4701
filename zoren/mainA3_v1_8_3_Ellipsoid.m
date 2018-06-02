%% Space Engineering 3 Assignment 3
% Author: Zoren Liu
% SID: 440353416
%

clear all; close all; clc;

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

% E_rot = Exx + Eyy + Ezz;

% e)
% Copper wire
l_Cu = 17.7; % cm
d_Cu = 1; % cm
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

A = [-a/2 -c/2 -b/2];
B = [a/2 -c/2 -b/2];
C = [-a/2 c/2 -b/2];
D = [-a/2 -c/2 b/2];
E = [-a/2 c/2 b/2];
F = [a/2 -c/2 b/2];
G = [a/2 c/2 -b/2];
H = [a/2 c/2 b/2];

origin = [0, 0, 0];

sim_type = 1;
axis_type = 1;

switch sim_type
    case 1
        P = [A;B;F;H;G;C;A;D;E;H;F;D;E;C;G;B];
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
        omega_x = -0.01;
        omega_y = 4;
        omega_z = 0.01;
    case 3 % About z
        omega_x = -0.01;
        omega_y = 0.01;
        omega_z = 4;
end

time = 60;
dt = 0.01;
tvec = (dt:dt:time)';

ang_x = 0;
ang_y = 0;
ang_z = 0;

counter = 0;

for n = 1:length(tvec)
    
    n
    
    omega_x_dot = (Iyy-Izz)*omega_y*omega_z/Ixx;
    omega_y_dot = (Izz-Ixx)*omega_z*omega_x/Iyy;
    omega_z_dot = (Ixx-Iyy)*omega_x*omega_y/Izz;
    
    omega_x = omega_x + omega_x_dot*dt;
    omega_y = omega_y + omega_y_dot*dt;
    omega_z = omega_z + omega_z_dot*dt;
    
    omegaDB(n, :) = [omega_x, omega_y, omega_z];
    
    ang_y = ang_y + omega_y*dt;
    ang_x = ang_x + omega_x*dt;
    ang_z = ang_z + omega_z*dt;
    
    angDB(n, :) = [ang_x, ang_y, ang_z];
    
    counter = counter + 1;
    
    E_rot(n, 1) = 1/2*Ixx*omega_x^2;
    E_rot(n, 2) = 1/2*Iyy*omega_y^2;
    E_rot(n, 3) = 1/2*Izz*omega_z^2;
    
    aElips1 = 2*sum(E_rot(n, :))/Ixx;
    bElips1 = 2*sum(E_rot(n, :))/Iyy;
    cElips1 = 2*sum(E_rot(n, :))/Izz;
    
    L(n, 1) = omega_x*Ixx;
    L(n, 2) = omega_y*Iyy;
    L(n, 3) = omega_z*Izz;
    
    Lsum(n) = sum(L(n,:));
    
    aElips2 = sum(L(n, :))^2/Ixx^2;
    bElips2 = sum(L(n, :))^2/Iyy^2;
    cElips2 = sum(L(n, :))^2/Izz^2;
    
    [x1, y1, z1] = ellipsoid(0,0,0,aElips1,bElips1,cElips1,30);
    [x2, y2, z2] = ellipsoid(0,0,0,aElips2,bElips2,cElips2,30);
    
    %     [h1, hel1] = imsurf(@(wx,wy,wz) wx.^2/aElips1+wy.^2/bElips1+wz.^2/cElips1, [-100,100,-100,100,-100,100]);
    %       hold on;axis equal;
    % [h2, hel2] = imsurf(@(wx,wy,wz) wx.^2/aElips2+wy.^2/bElips2+wz.^2/cElips2, [-100,100,-100,100,-100,100]);
    %       h2.FaceAlpha = 0.6;
    %       h = intercurve(hel1, hel2);
    %       h.Color = 'r';
    %       camlight;
    
    
    %     [intMtr, intSurf] = SurfaceIntersection([x1, y1, z1], [x2, y2, z2]);
    %     pause
    
%     figure(1)
%     if counter == 10
%         counter = 0;
% 
%         switch sim_type
%             case 1
%                 
%                 dcmMtr = dcm([ang_y; ang_x; ang_z]);
%                 P_plot = P*dcmMtr;
%                 subplot(1,2,1)
%                 plot3(P_plot(:,1),P_plot(:,2),P_plot(:,3),'k', 'LineWidth', 2), hold on % original cube
%                 plot3(P_plot(1:2,1),P_plot(1:2,2),P_plot(1:2,3),'b', 'LineWidth', 2)
%                 plot3(P_plot(8:9,1),P_plot(8:9,2),P_plot(8:9,3),'r', 'LineWidth', 2)
%                 plot3(P_plot(2:3,1),P_plot(2:3,2),P_plot(2:3,3),'g', 'LineWidth', 2)
%                 plot3(0,0,0,'or')
%                 axis equal; grid on; axis([-0.20 0.20 -0.20 0.20 -0.20 0.20]); rotate3d on; hold off
%                 drawnow
%                 
%             case 2
%                 
%                 subplot(1,2,1)
%                 rotate(P, [0 1 0], rad2deg(omega_x*dt*10), origin);
%                 rotate(P, [1 0 0], rad2deg(omega_y*dt*10), origin);
%                 rotate(P, [0 0 1], rad2deg(omega_z*dt*10), origin);
%                 axis equal; grid on; axis([-0.20 0.20 -0.20 0.20 -0.20 0.20]); rotate3d on; hold off
%                 drawnow
%                 
%         end
%         
%         subplot(1,2,2)
%         plot(tvec(1:n), euler(1, 1:n), tvec(1:n), euler(2, 1:n), tvec(1:n), euler(3, 1:n),...
%             tvec(1:n), Lsum(1:n),'LineWidth', 1.5)
%         xlabel('Time (s)'); ylabel('\omega (rad/s)');
%         legend('\omega_x', '\omega_y', '\omega_z')
%         grid on
%         
%         
%         figure(2)
%         surf(x1, y1, z1, 'edgecolor', 'none','FaceAlpha',0.3,'FaceColor',[0.9100 0.4100 0.1700])
%         axis([-250 250 -250 250 -250 250]); axis equal; hold on
%         surf(x2, y2, z2, 'edgecolor', 'none','FaceAlpha',0.3,'FaceColor',[0 0 1])
%         drawnow
%         axis([-250 250 -250 250 -250 250]); axis equal; view(-45, 45);  grid on;hold off;
%         
%     end
    
    
    
    
end

figure(2)
plot(tvec, omegaDB(:, 1), tvec, omegaDB(:, 2), tvec, omegaDB(:, 3),...
    'LineWidth', 1.5)
xlabel('Time (s)'); ylabel('\omega (rad/s)');
legend('\omega_x', '\omega_y', '\omega_z')
grid on
figure(3)
plot(tvec, angDB(:, 1), tvec, angDB(:, 2), tvec, angDB(:, 3),...
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



