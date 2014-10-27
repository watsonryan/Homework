function deltacr = gendelcr3(lat,vel_L,height,deltat,DCMbn,DCMel,earthflg)
%GENDELCR3   Function to generate the component of delta-theta
%            associated with craft-rate at a single point on a path
%            profile generated with velocity in ENU coordinates
%
%	deltacr = gendelcr3(lat,vel_L,height,deltat,DCMbn,DCMel,earthflg)
%
%   INPUTS
%       lat = latitude of the current point (radians)
%       vel_L = velocity in ENU coordinates (meters/sec)
%       height = vehicle height above the reference 
%                     ellipsoid (meters) for each path segment
%       deltat = time interval between data points (seconds)
%       DCMbn = body-to-nav direction cosine matrix
%       DCMel = earth-frame to local-level-frame direction 
%               cosine matrix
%       earthflg = earth shape flag
%                  0 = spherical earth; 1 = WGS-84 ellipsoidal earth
%                  (see CRAFRATE)
%
%   OUTPUTS
%       deltacr(1,1:3) = the three delta-theta
%                        components (in the body frame) of craft rate
%

%	M. & S. Braasch 6-98; 3-05; 11-08
%	Copyright (c) 1997-2008 by GPSoft LLC
%	All Rights Reserved.
%
%   REVISION HISTORY
%
%   November 2008:  created the function to work on a single data point

if nargin<7,error('insufficient number of input arguments'),end

C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

   vx = vel_L(1);
   vy = vel_L(2);
   rho = crafrate(lat,vx,vy,height,DCMel,earthflg,0);
   rho_b = (DCMbn')*(C*rho);     % convert from local-level to body coordinates
   deltacr(1,1:3) = ( rho_b*deltat )';   
   

