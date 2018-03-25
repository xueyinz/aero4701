function [m, xhat, xUKF, qIntO] = adcs(sensor,dt,sw,xhat,bo)

    % Extract data structures
    cont = sw.cont;
    esti = sw.esti;
    
    % Get sun sensor vector
    sb = sunSensor(sensor.ss);
    
    % Get the acceleromater data
    [b1, w1] = IMU(sensor.IMU1);
    [b2, ~] = IMU(sensor.IMU2);
    
    % Create the measurement array
    y = [w1; b1; b2; sb];
    
    % Average the magnetic field readings for use in the controller
    b = (b1 + b2)/2;        % If time, look into estimating this too
    
    % Initialise persistent variables the first time around - this
    % simulates the memory of the ARDUINO
    persistent xhatUKF Pxx  mode bm1 qInt bdotm1 GPS_last
    if isempty(xhatUKF)
        xhatUKF = [0; 0; 0; esti.bias0];
        Pxx = sw.esti.P0;
        mode = sw.cont.mode;
        bm1 = b;
        bdotm1 = zeros(3,1);
        qInt = zeros(3,1);
        GPS_last = sensor.GPS;
        fprintf('Initialisation correct\n\n');
    end
    
    % NOTE: xhat should appear in the list above, but has been extracted
    % for debugging
    if all(isnan(xhat))
        xhat = [w1 + esti.bias0; esti.q0];
    end
    
    % Calculate the ECI to orbit DCM
    [Rio, ~, wo] = orbitProp(sensor.GPS, GPS_last, dt);
    GPS_last = sensor.GPS;
    esti.wo = wo;
    cont.wo = wo;
    
    % Calculate what the magnetic field and sun vector should be
    % Note: for bo we use the c library here, so we just use Matlab's
    % wrldmagm in this simulation
    % https://www.ngdc.noaa.gov/geomag/WMM/thirdpartycontributions.shtml
    [~, so] = vecModel(Rio, sensor.GPS);

    % Estimate the state matrix - note the b needs to be changed for an
    % estimation of what it should be using a simple dipole model
    [xhatUKF, Pxx, xhat] = UKF(esti, xhatUKF, Pxx, y, xhat, bo, so, dt);
    
    % For output
    xUKF = xhatUKF;
    
    % Integrate the quaternions
    qInt = qInt + xhat(4:6)*dt;
    
    % Mode selector
    [mode,qInt] = modeSelect(mode,xhat,cont,qInt);
    qIntO = qInt;
    
    % Calcualte the control
    [m, bdot] = controlLQR(xhat,qInt,mode,cont,b,bm1,bdotm1);
    bdotm1 = bdot;
    bm1 = b;
    
    % Saturate outputs
    m(m > cont.mmax) = cont.mmax(m > cont.mmax);
    m(m < -cont.mmax) = cont.mmax(m < -cont.mmax);