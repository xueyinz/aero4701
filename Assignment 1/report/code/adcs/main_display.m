%% The MuirKat - ADCS Subsystem
% Simulation for system validation
% Mathew Gardiner
clc; clear; clear functions; clear all; %#ok<*CLALL,CLFUNC>

% NOTE: sw structure contains variables accesible by the Arduino code

% Add folders
addpath('Dynamics');
addpath('Other');
addpath('Software');
addpath('Hardware');

%% Sim Settings

% Time
orbits = 6;
dt = 10;

% Random seed - for report results
rng(65);

% Inital state
eul0 = randi([-180,180],3,1)*pi/180;
pqr0 = [10; 10; 10]*pi/180;

% Initial time
sat.dvec0 = [2018 6 30 12 0 0];

%% Estimator Parameters

% Initial conditions
sw.esti.q0 = E2Q([0; 0; 0]);            % Initial quaternion estimate
sw.esti.bias0 = [0; 0; 0]*pi/180;       % Initial bias estimate

% Initial covariance
sw.esti.P0 = [(50*pi/180)^2*eye(3) zeros(3);
    zeros(3) (0.2*pi/180)^2*eye(3)];

% Other UKF parameters
sw.esti.a = 1;                          % Scaling parameter
sw.esti.lambda = 1;                     % Controls number of sigma points

%% Sensor parameters

% Noise
sensor.var.v = (0.1*pi/180)^2;          % Noise on gyro
sensor.var.u = (1.75e-6)^2;             % Noise on gyro drift
sensor.var.b = (0.5e-6)^2;              % Noise on magnetometer
sensor.var.a = (8e-3/9.81)^2;           % Noise on accelerometer
sensor.var.X = 25^2;                    % Noise on GPS position
sensor.var.V = 1^2;                     % Noise on GPS velocity

% Number of samples per time step
sensor.samples = 100;

% Gyroscope bias
sensor.dbias = (2*rand(3,1)-1)*pi/180/3600;     % Bias drift, deg/hour
sensor.bias0_true = [0.1; 0.1; 0.1]*pi/180;     % Bias noise

% Estimator noise - calculated
sw.esti.var.v = sensor.var.v/sensor.samples;    % Noise on gyro
sw.esti.var.u = sensor.var.u/sensor.samples;    % Noise on gyro drift
sw.esti.var.b = sensor.var.b/sensor.samples;    % Noise on magnetometer
sw.esti.var.a = sensor.var.a/sensor.samples;    % Noise on accelerometer

%% True Satellite Parameters

% Physical parameters
sat.I = [354277.7 10941.55 -3809.14;        % Moment of inertia
    10941.55 375364.21 -1870.74;
    -3809.14 -1870.74 373076]/(1000*1000*1000);

% Rotate moment of inertia matrix from M to B frame
sat.I = (DCM(pi,1)*DCM(pi/2,2))*sat.I*(DCM(pi,1)*DCM(pi/2,2))';

% Orbit parameters
orb.altk = 400;                             % Altitude (km)
orb.i = 51.6413;                            % Inclination (deg)
orb.RA = 360*rand(1);                       % Right Ascension (deg)
orb.AP = 360*rand(1);                       % Argument of the perigee (deg)

% Drag parameters
drag.on = true;                             % Drag switch
drag.Cd = 2.2;                              % Coefficient of drag - cube
drag.CdL = 1.17;                            % Coefficient of drag - cyl
drag.rho = 3.89e-12;                        % Air density
drag.faces = [0.1^2; 0.1*0.115; 0.1*0.115]; % Face areas
langmuirL = 45;                             % Probe protrusion, mm
drag.langVec = [0; 0; 1];                   % Langmuir probe direction
drag.langmuirS = langmuirL*2.33e-3/1000;    % Probe area
X = 49.89; Y = -46.98; Z = 56.11;           % Location of CG (solidworks)

% Locations of each face
faces = [115/2, 0, 0; -115/2, 0, 0; 0, 50, 0;
    0, -50, 0; 0, 0, 50; 0, 0, -50];

% Calculate the centre of pressure of the faces and the probes, converting
% from the SolidWorks axis to body axis
drag.rcp = (faces - repmat([-(Z-115/2), -(Y+50), -(X-50)],6,1))';
drag.rcpL = ([115/2 28 50+langmuirL/2; 115/2 -28 -50-langmuirL/2] - ...
    repmat([-(Z-115/2), -(Y+50), -(X-50)],2,1))';
drag.rcp = drag.rcp/1000; drag.rcpL = drag.rcpL/1000;

%% Controller Software Settings

% Mode
sw.cont.mode = 'bdot';     % Starting mode, 'bdot' or 'point' or 'none'

% Bdot gain
sw.cont.bdotK = 3.4e-6;
cont.bdotwc = 0.7;          % Filtering frequency

% Conditions to switch modes
sw.cont.endBdot = 0.1*pi/180;               % End detumbling, rad/s
sw.cont.startBdot = 2*sw.cont.endBdot;      % Start detumbling, rad/s

% Gain matrix
cont.qq = 1/0.05^2;                 % Quaternion weight
cont.qw = 1/(0.04*pi/180)^2;        % Rotation weight
cont.rmax = 1/0.05^2;               % Control weight
sw.cont.Ki = 1e-4;                  % Integrator gain

% Time step for initial averaging of B matrix
cont.dt = 10;
cont.Norb = 1;

% Maximum control ouputs and coil parameters
coil.N = [150; 250; 100];                           % Coil turns
coil.A = [(22/2000)^2*pi; (22/2000)^2*pi; 0.08^2];  % Coil areas
coil.NA = coil.N.*coil.A;                           % Area and turns
coil.R = [5.277969906; 8.79661651; 21.7216];        % Coil resistance
sw.cont.Imax = [0.2; 0.2; 0.2];                     % Maximum current
sw.cont.mmax = sw.cont.Imax./coil.NA;               % Maximum moment
sw.cont.coilNA = coil.NA;

%% Further Orbit Calculations

% Orbit initialisations
orb.alt = orb.altk*1000;                        % Altitude
orb.mu = 3.986004418e14;                        % Gravitational parameter
orb.Re = 6371e3;                                % Earth radius
orb.i = orb.i*pi/180;                           % Inclination
orb.T = 2*pi*sqrt((orb.alt+orb.Re)^3/orb.mu);   % Period
orb.V = sqrt(orb.mu/(orb.Re+orb.alt));          % Velocity
orb.wo = sqrt(orb.mu/(orb.Re+orb.alt)^3);       % Angular velocity
orb.dvec0 = sat.dvec0;                          % Date vector

%% Calculate the controller gain

% Estimated moment of inertia matrix
cont.I = sat.I;

% Controller weights
cont.Q = diag([ones(1,3)*cont.qw,ones(1,3)*cont.qq]);
cont.R = eye(3)*cont.rmax;

% Useful parameters
Ix = cont.I(1,1); Iy = cont.I(2,2); Iz = cont.I(3,3);
sigx = (Iy-Iz)/Ix; sigy = (Iz - Ix)/Iy; sigz = (Ix - Iy)/Iz;

% Calculate the A matrix
A = zeros(6);
A(1,3) = -orb.wo*sigx; A(3,1) = -orb.wo*sigz;
A(1,4) = -6*orb.wo^2*sigx; A(2,5) = 6*orb.wo^2*sigy;
A(4:6,1:3) = 0.5*eye(3);
A(4,6) = -orb.wo; A(6,4) = orb.wo;
cont.A = A;

% Calculate or load the B matrix
% cont = ABcalc(cont,orb);
% B = cont.B;
% save('Bmat','B');
data = load('Bmat');
cont.B = data.B;

% Find the constant gain matrix
sw.cont.K = lqr(cont.A,cont.B,cont.Q,cont.R);

% Stability check
fprintf('Eigenvalues of closed loop system:\n\n')
disp(vpa(eig(cont.A-cont.B*sw.cont.K),3));

% Set up the filter
cont.H = tf([cont.bdotwc 0], [1 cont.bdotwc]);
cont.Hd = c2d(cont.H,dt,'tustin');
[sw.cont.tf.num, sw.cont.tf.den] = tfdata(cont.Hd);

%% Initialise the simulation

% Rotation matrix function handle
R = @(q) [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2, ...
    2*(q(1)*q(2) + q(3)*q(4)), 2*(q(1)*q(3) - q(2)*q(4));
    2*(q(1)*q(2) - q(3)*q(4)), -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,...
    2*(q(2)*q(3) + q(1)*q(4));
    2*(q(1)*q(3) + q(2)*q(4)), 2*(q(2)*q(3) - q(1)*q(4)), ...
    -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2];

% Things for the other structures
sat.i = orb.i; sat.T = orb.T; sat.wo = orb.wo; sat.drag = drag;
cont.wo = orb.wo; esti.wo = orb.wo;

% Time steps
T = orb.T*orbits;
tvals = 0:dt:T;
N = numel(tvals);

% Preallocate arrays
m = zeros(3,N); lla = m; b = m; wo = m; qInt = m; biasVals = m;
Tcvals = m; Tdvals = m;
xLin = zeros(6,N); xUKF = xLin;
xAng = zeros(7,N); xhat = xAng;

% Initial conditions
xLin(:,1) = orbitInit(orb.altk,orb.altk,orb.i,orb.RA,orb.AP);
quats = E2Q(eul0);
wo(:,1) = pqr0;
pqrb0 = pqr0 + R(quats)*[0; orb.wo; 0];
xAng(:,1) = [quats; pqrb0];
X = [xLin; xAng];
[b(:,1), bo, lla(:,1)] = bBody(X(:,1),0,orb);
xhat(:,1) = NaN(7,1);
bias1 = repmat(sensor.bias0_true,1,sensor.samples);
bias2 = zeros(3,sensor.samples);

%% Loop through time
fprintf('Beginning Simulation\n\n'); tic;
for idx = 2:N
    
    % Time
    t = tvals(idx);
    
    % Body dynamics and integrate
    X(:,idx) = RK4(@satDyn,X(:,idx-1),m(:,idx-1),dt,b(:,idx-1),sat);
    
    % Normalise the Quaternions
    X(7:10,idx) = X(7:10,idx)./(norm(X(7:10,idx)));
    
    % Calculate b in the control frame
    [b(:,idx), lla(:,idx), wo(:,idx), bo] = bBody(X(:,idx),t,orb);
    
    % Sun-sensor outputs (including noise)
    sensorOut.ss = sunSensor_phys();
    
    % IMU ouputs (including noise)
    [sensorOut.IMU1, bias1] = IMU_phys(sensor,bias1,X(:,idx),b(:,idx),dt);
    [sensorOut.IMU2, bias2] = IMU_phys(sensor,bias2,X(:,idx),b(:,idx),dt);
    biasVals(:,idx) = mean(bias1,2);
    
    % GPS outputs
    sensorOut.GPS = GPS_phys(X(:,idx),t,orb,sensor);
    
    % Software loop
    [m(:,idx), xhat(:,idx),xUKF(:,idx),qInt(:,idx)] = ...
        adcs(sensorOut,dt,sw,xhat(:,idx-1),bo);
    
end

%% Helper functions
function [b, lla, wo, bo] = bBody(X,t,orb)
% Returns b in control frame

% Rotation matrix function handle
R = @(q) [q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2, ...
    2*(q(1)*q(2) + q(3)*q(4)), 2*(q(1)*q(3) - q(2)*q(4));
    2*(q(1)*q(2) - q(3)*q(4)), -q(1)^2 + q(2)^2 - q(3)^2 + q(4)^2,...
    2*(q(2)*q(3) + q(1)*q(4));
    2*(q(1)*q(3) + q(2)*q(4)), 2*(q(2)*q(3) - q(1)*q(4)), ...
    -q(1)^2 - q(2)^2 + q(3)^2 + q(4)^2];

% Latitude and longitude
lla = eci2lla(X(1:3)',dTcalc(orb.dvec0,t));

% Calcualte the b vector in ECI
bned = wrldmagm(lla(3)/1000, lla(1), lla(2), 2018.5)*10^-9;
beci = (dcmeci2ecef('IAU-2000/2006',dTcalc(orb.dvec0,t))\...
    (dcmecef2ned(lla(1), lla(2))\bned));

% Find the quaternion in ECI
qeci = [sin(orb.i/2)*cos((orb.RA*pi/180 - orb.wo*mod(t,orb.T))/2);
    sin(orb.i/2)*sin((orb.RA*pi/180 - orb.wo*mod(t,orb.T))/2);
    cos(orb.i/2)*sin((orb.RA*pi/180 - orb.wo*mod(t,orb.T))/2);
    -cos(orb.i/2)*cos((orb.RA*pi/180 - orb.wo*mod(t,orb.T))/2)];

% Change frames
Ro2eci = R(qeci); Ro2c = R(X(7:10));
bo = Ro2eci\beci;
bo = [bo(2); -bo(3); -bo(1)];
b = Ro2c*bo;

% Find rotation in the orbit frame
wo = X(11:13) - R(X(7:10))*[0; orb.wo; 0];

end