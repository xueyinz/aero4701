clear;
clc; 
close all; 

% This initializes all constants as global. 
constants; 

% Solver settings
dt = 0.01; 
tMax = 10; 
noSamples = ceil(tMax/dt); 

% Always initialize space for the solution. 
x = NaN(10, noSamples); 

% Initial quaternion 
x(1:4, 1) = [1; 0; 0; 0];

% Initial angular velocity
x(5:7, 1) = [pi/2; 0; 0]; 

% Initial reaction wheel velocity
x(8:10, 1) = [0; 0; 0]; 

% functionHandle = @(arg1, arg2) ( arg1 + arg2 ); 

% Initialize the function handle. 
solveHandle = @(x)( f(x, controller(x)) ); 

% Initialize waitbar
h = waitbar(0, 'Solving the problem...'); 
for i = 2:noSamples
    
    % Solve.
    x( :, i) = RK4(solveHandle, x(:,i-1), dt); 
    % Update waitbar
    waitbar(i/noSamples, h); 
    
end
% Delete waitbar
delete(h); 

time = linspace(0, tMax, noSamples); 
skip = 5; 
% Animation
for i = 1:skip:noSamples
    % trajectory3 is not my code - hacked it to use quaternions. 
    trajectory3(0, 0, 0, 0, 0, 0, x(1:4, i)', 1, 0.1, 'shuttle');
    % drawnow() flushes the graphics buffer. Updates the graphics at each
    % iteration
    drawnow(); 
    
end

figure(2)
subplot(3,1,1);
plot(time, x(5, :));
ylabel('w_s_x');
title('Spacecraft Angular Velocity'); 
subplot(3,1,2);
plot(time, x(6, :));
ylabel('w_s_y'); 
subplot(3,1,3); 
plot(time, x(7, :)); 
ylabel('w_s_z'); 
xlabel('time (s)'); 


figure(3)
subplot(3,1,1);
plot(time, x(8, :));
ylabel('w_w_x');
title('Reaction Wheel Angular Velocity'); 
subplot(3,1,2); 
plot(time, x(9, :));
ylabel('w_w_y'); 
subplot(3,1,3);
plot(time, x(10, :));
ylabel('w_w_z'); 
xlabel('time (s)'); 
