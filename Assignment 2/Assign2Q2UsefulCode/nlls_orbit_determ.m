function params = nlls_orbit_determ(obs,GS_ECEF,init_posvel_guess)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NLLS_ORBIT_DETERM - Uses Non-linear Least Squares technique to estimate
% the orbital parameters given a series of ground station tracking
% obserations.
%
% Inputs: obs - N x 4 vector of [time, range, azimuth, elevation] measurements
%         from a single ground station
%         GS_ECEF - ECEF position of the ground tracking station
%         init_posvel_guess - Initial parameter estimate 
%
% Outputs: params - Orbital parameters: [a,e,i,RAAN,AoP,Mo (at epoch)]
%

% NOTE: Elements of this function are missing! you will need to fill in
% these gaps to get the code working!

% This version of the code is designed for a single ground station; you
% will need to make several adjustments for multiple ground stations 

% See Week 6 Slides 18 to 30 for Details on non-linear least squares for
% orbit determination from ground station tracking measurements

% FILL IN HERE
time_last_vernal_equinox = ....

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialise H matrix and error vector

num_obs = size(obs,1); % Calculate the total number of ground station observations
H = zeros(3*num_obs,6); % Initialise the total measurement Jacobian matrix
delta_y = zeros(3*num_obs,1); % Initialise the residual vector

% Use initial guess at orbital parameters to get initial position and
% velocity 
X = init_posvel_guess;
delta_x = 10e9;

% Main iteration loop to converge on initial satellite position and
% velocity
max_iter = 20;
tol = 1e-6;
iter = 0;
while (norm(delta_x) > tol)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Re-linearise the system, build the observation model and Jacobian
    
    % Loop through each observation, build observation matrix
    for i = 1:num_obs
        
        % Use current value of predicted initial position to calculate
        % position at time of observation
        current_obs_time = obs(i,1);
        [pos, vel] = universal_conic_section_orbit(current_obs_time - time_last_vernal_equinox, X(1:3), X(4:6));
        
        % Get the expected measurement (range, azimuth, elevation)
        % FILL IN HERE: expected measurement as a function of ground
        % station location, and satellite location
        expected_obs = ...
        
        % Build H matrix
        dxdxo = universal_conic_section_jacobian(pos, vel, X(1:3), X(4:6), current_obs_time - time_last_vernal_equinox);
        dhdx = rangeazielev_obs_jacob(GS_ECEF, pos, current_obs_time - time_last_vernal_equinox);
        H((3*i - 2):3*i,:) = dhdx*dxdxo; % Builds the rows of the Jacobian matrix corresponding to this measurement
        
        % Build the residual vector for this measurement
        actual_obs = obs(i,2:4);
        delta_y((3*i - 2):3*i,1) = (actual_obs - expected_obs)';
        
    end
    
    % Use non-linear least squares to estimate error in x
    delta_x = inv(H'*H)*H'*delta_y;
    
    % Alternative formulation (much better, type "help mldivide" for
    % details)
    % delta_x = (H'*H)\H'*delta_y;
    
    X = X + delta_x;
    
    iter = iter + 1;
    
    if (iter >= max_iter)
        disp('Failed to Converge !!')
        break;
    end
    
end

% FILL IN HERE: Transform iterated estimate of initial position and
% velocity into orbital parameters (See Notes Week 6 Slide 30)

% .......

% end of nlls_orbit_determ






