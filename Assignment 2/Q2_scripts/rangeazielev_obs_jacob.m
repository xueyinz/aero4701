function dhdx = rangeazielev_obs_jacob(GS_ECEF, GS_LLH, sat_ECI, time_svep)

% RANGEAZIELEV_OBS_JACOB - computes the jacobian matrix for the function to
% convert a sat and GS position and time into a range, azimuth and
% elevation observation
%
% Inputs: GS_ECEF - The ground station location in ECEF coordinates
%         sat_ECI - The ECI position of the satellite at the observation
%         time_svep - Time since the last vernal equinox passage
%
% Outputs: dhdx - Jacobian matrix of observation function w.r.t position
% and velocity at the observation time (3x6 matrix)
%

% NOTE: Elements of this function are missing! you will need to fill in
% these gaps to get the code working!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILL IN HERE: Functions for computing "Ceci2ecef" (ECI to ECEF DCM),
% function for computing "Clgcv2ecef" (LGCV to 

Ceci2ecef = C_eci2ecef(time_svep);
    
Clgcv2ecef = C_lgcv2ecef(GS_LLH);

% compute NED cartesian range vector
satrelNED = Clgcv2ecef'*(Ceci2ecef*sat_ECI - GS_ECEF);
dx = satrelNED(1);
dy = satrelNED(2);
dz = satrelNED(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jacobian for the transformation from cartesian (x,y,z) to polar (r,a,e)
% (See Week 6 Notes)
%
% Note: You may wish to have different forms of observations (i.e. not just
% z = [Range, Azimuth and Elevation] but perhaps Angles-Only (i.e.
% z = [Azimuth, Elevation]) or Range only (i.e. z = [Range]). You can vary
% the jacobian format in this function if desired.
%

rxy2 = dx^2+dy^2;
rxy  = sqrt(rxy2);
r2 = dx^2+dy^2+dz^2;
r  = sqrt(r2);

Hc2p(1,1) = dx/r;       		Hc2p(1,2) = dy/r;       		Hc2p(1,3) = dz/r;
Hc2p(2,1) = -dy/rxy2;   		Hc2p(2,2) = dx/rxy2;    		Hc2p(2,3) = 0;
Hc2p(3,1) = -dx*dz/(r2*rxy);	Hc2p(3,2) = -dy*dz/(r2*rxy); 	Hc2p(3,3) = rxy/r2;

% dhdx
dhdx = [Hc2p*Clgcv2ecef'*Ceci2ecef, zeros(3)];

% end of rangeazielev_obs_jacob
