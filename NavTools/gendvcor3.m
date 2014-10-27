function dvcor = gendvcor3(lat,vel_L,height,deltat,DCMbn,DCMel,earthflg)
%GENDVCOR3      Function to generate the component of delta-V
%               associated with Coriolis and gravity.  This version
%               provides a value for a single sample point from
%               a velocity profile which has been 
%               generated in ENU coordinates.  North-pointing mechanization
%               assumed.
%       
%	dvcor = gendvcor3(lat,vel_L,height,time,DCMbn,DCMel,earthflg)
%
%   INPUTS
%       lat = latitude (radians) of the sample point from the flight path
%       vel_L = velocity in ENU (m/s)
%       height = vehicle ellipsoidal height (meters)
%       time = time interval between sample points (seconds)
%       DCMbn = body-to-nav profile of direction cosine matrix
%       DCMel = earth-frame to local-level frame direction cosine matrix
%       earthflg = earth shape flag
%                  0 = spherical earth; 1 = WGS-84 ellipsoid (see CRAFRATE)
%
%   OUTPUTS
%       dvcor = Coriolis and gravity components of delta-V
%               in body frame (nose-rt.wing-down)
%

%	M. & S. Braasch 6-98; 3-05; 11-08
%	Copyright (c) 1998-2008 by GPSoft
%	All Rights Reserved.
%
%   REVISION HISTORY
%
%   November 2008:  function created to operate on a single sample point

if nargin<7,error('insufficient number of input arguments'),end
if (earthflg ~= 0) && (earthflg ~= 1), error('EARTHFLG not specified correctly'),end

C = [0 1 0; 1 0 0; 0 0 -1];      % Conversion between ENU and NED
vertmech = 0;

   vx = vel_L(1);
   vy = vel_L(2);
  
   omega_el_L = crafrate(lat,vx,vy,height,DCMel,earthflg,vertmech);
   omega_en_n = C*omega_el_L;
   
   vel_cor_b = coriolis(vel_L,omega_en_n,DCMbn,DCMel,1,deltat);
  
   g_tru = gravity(lat,height);
   g_vect_n = [0 0 -g_tru]';
   vel_g_n = g_vect_n*deltat;

   vel_g_b = (DCMbn')*vel_g_n;
  
   dvcor(1,1:3) = (-vel_cor_b + vel_g_b)';


