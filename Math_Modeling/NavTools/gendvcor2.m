function dvcor = gendvcor2(lat_prof,vel_prof_L,height_prof,time,...
                          DCMnb_prof,DCMel_prof,earthflg)
%GENDVCOR2      Function to generate the component of delta-V
%               associated with Coriolis and gravity.  This version
%               operates on a velocity profile which has been 
%               generated in ENU coordinates.  North-pointing mechanization
%               assumed.
%       
%	dvcor = gendvcor2(lat_prof,vel_prof_L,height_prof,time,...
%                          DCMnb_prof,DCMel_prof,earthflg)
%
%   INPUTS
%       lat_prof = profile of latitude (radians) of the flight path
%       vel_prof_L = velocity profile in ENU (m/s)
%       height_prof = vehicle height profile over the flight path (meters)
%       time = sequential time vector over the flight path (seconds)
%       DCMnb_prof = profile of direction cosine elements over time
%                    relating nav-frame to body-frame
%          DCMnb_prof(i,1:9) = elements of the i-th direction cosine matrix
%                            (DCM) for vehicle attitude (navigation-to-
%                            body); 1 = DCM(1,1),
%                            2 = DCM(1,2), 3 = DCM(1,3),
%                            4 = DCM(2,1), et cetera
%       DCMel_prof = profile of direction cosine elements over time
%                    relating earth-frame to local-level-frame
%          DCMel_prof(i,1:9) = elements of the i-th direction cosine matrix
%                            (DCM) for vehicle position; 1 = DCM(1,1),
%                            2 = DCM(1,2), 3 = DCM(1,3),
%                            4 = DCM(2,1), et cetera
%       earthflg = earth shape flag
%                  0 = spherical earth; 1 = WGS-84 ellipsoid (see CRAFRATE)
%
%   OUTPUTS
%       dvcor = profile of delta-V (Coriolis and gravity only)
%               in body frame (nose-rt.wing-down)
%

%	M. & S. Braasch 6-98; 3-05
%	Copyright (c) 1998-2005 by GPSoft
%	All Rights Reserved.
%

if nargin<7,error('insufficient number of input arguments'),end
if (earthflg ~= 0) && (earthflg ~= 1), error('EARTHFLG not specified correctly'),end

C = [0 1 0; 1 0 0; 0 0 -1];      % Conversion between ENU and NED
vertmech = 0;

for i = 2:length(lat_prof),

   vx1 = vel_prof_L(i-1,1);
   vy1 = vel_prof_L(i-1,2);
   vz1 = vel_prof_L(i-1,3);

   vx2 = vel_prof_L(i,1);
   vy2 = vel_prof_L(i,2);
   vz2 = vel_prof_L(i,3);
   td12 = time(i) - time(i-1);
   tdin = 0.5*td12;
   lat_in = interpol(lat_prof(i-1),lat_prof(i),td12,tdin);
   vx_in = interpol(vx1,vx2,td12,tdin);
   vy_in = interpol(vy1,vy2,td12,tdin);
   vz_in = interpol(vz1,vz2,td12,tdin);
   vel_in = [vx_in vy_in vz_in];
   height_in = interpol(height_prof(i-1),height_prof(i),td12,tdin);

   DCMnb2=[DCMnb_prof(i,1:3); DCMnb_prof(i,4:6); DCMnb_prof(i,7:9)];
   DCMnb1=[DCMnb_prof(i-1,1:3); DCMnb_prof(i-1,4:6); DCMnb_prof(i-1,7:9)];
   DCMnbavg = 0.5*( DCMnb2 + DCMnb1 );
   DCMbnavg = DCMnbavg';
   DCMel2=[DCMel_prof(i,1:3); DCMel_prof(i,4:6); DCMel_prof(i,7:9)];
   DCMel1=[DCMel_prof(i-1,1:3); DCMel_prof(i-1,4:6); DCMel_prof(i-1,7:9)];
   DCMelavg = 0.5*( DCMel2 + DCMel1 );
  
   omega_el_L = crafrate(lat_in,vx_in,vy_in,height_in,DCMelavg,earthflg,vertmech);
   omega_en_n = C*omega_el_L;
   
   vel_cor_b = coriolis(vel_in,omega_en_n,DCMbnavg,DCMelavg,1,td12);
  
   g_tru = gravity(lat_in,height_in);
   g_vect_n = [0 0 -g_tru]';
   vel_g_n = g_vect_n*td12;

   vel_g_b = DCMnbavg*vel_g_n;
  
   dvcor(i-1,1:3) = (-vel_cor_b + vel_g_b)';
end

