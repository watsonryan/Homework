function deltacr = gendelcr2(lat_prof,vel_prof_L,height_prof,...
    time_prof,DCMnb_prof,DCMel_prof,earthflg)
%GENDELCR2   Function to generate the component of delta-theta
%            associated with craft-rate for a path
%            profile generated with velocity in ENU coordinates
%
%	deltacr = gendelcr2(lat_prof,vel_prof_L,height_prof,...
%                       time_prof,DCMnb_prof,DCMel_prof,earthflg)
%
%   INPUTS
%       lat_prof = vector of latitudes for each waypoint in the path
%       (radians)
%       vel_prof_L = velocity profile in ENU coordinates (meters/sec)
%       height_prof = vector of vehicle height above the reference 
%                     ellipsoid (meters) for each path segment
%       time_prof = time vector (seconds)
%       DCMnb_prof = profile of nav-frame (NED) to body-frame direction 
%                    cosine elements over time
%                    DCMnb_prof(i,1:9) = elements of the i-th direction cosine matrix
%                                   (DCM) for vehicle attitude; 1 = DCM(1,1),
%                                   2 = DCM(1,2), 3 = DCM(1,3),
%                                   4 = DCM(2,1), et cetera
%       DCMel_prof = profile of earth-frame to local-level-frame direction 
%                    cosine elements over time
%                    DCMel_prof(i,1:9) = elements of the i-th direction cosine matrix
%                                   (DCM) for vehicle position; 1 = DCM(1,1),
%                                   2 = DCM(1,2), 3 = DCM(1,3),
%                                   4 = DCM(2,1), et cetera
%       earthflg = earth shape flag
%                  0 = spherical earth; 1 = WGS-84 ellipsoidal earth
%                  (see CRAFRATE)
%
%   OUTPUTS
%       deltacr(i,1:3) = for the i-th path segment, the three delta-theta
%                        components (in the body frame) of craft rate
%

%	M. & S. Braasch 6-98; 3-05
%	Copyright (c) 1997-2005 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<7,error('insufficient number of input arguments'),end

%%nmph2mps = 1.6878*0.3048;
C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

for i = 1:( max(size(lat_prof)) - 1 ),
   %%velmps = totvel_prof(i)*nmph2mps;   %Convert from nautical-miles-per-hour to meters-per-second
   DCMnb1=[DCMnb_prof(i,1:3); DCMnb_prof(i,4:6); DCMnb_prof(i,7:9)];
   DCMnb2=[DCMnb_prof(i+1,1:3); DCMnb_prof(i+1,4:6); DCMnb_prof(i+1,7:9)];
   DCMnbavg = 0.5*( DCMnb1 + DCMnb2 );
   DCMel1=[DCMel_prof(i,1:3); DCMel_prof(i,4:6); DCMel_prof(i,7:9)];
   DCMel2=[DCMel_prof(i+1,1:3); DCMel_prof(i+1,4:6); DCMel_prof(i+1,7:9)];
   DCMelavg = 0.5*( DCMel1 + DCMel2 );
   %vx1 = velmps*sin(tc_prof(i));
   %vx2 = velmps*sin(tc_prof(i+1));
   vx1 = vel_prof_L(i,1);
   vx2 = vel_prof_L(i+1,1);
   vxavg = 0.5*(vx1 + vx2);
   %vy1 = velmps*cos(tc_prof(i));
   %vy2 = velmps*cos(tc_prof(i+1));
   vy1 = vel_prof_L(i,2);
   vy2 = vel_prof_L(i+1,2);
   vyavg = 0.5*(vy1 + vy2);
   latavg = 0.5*(lat_prof(i) + lat_prof(i+1));
   heightavg = 0.5*(height_prof(i) + height_prof(i+1));
   rho = crafrate(latavg,vxavg,vyavg,heightavg,DCMelavg,earthflg,0);
   rho_b = DCMnbavg*(C*rho);     % convert from local-level to body coordinates
   deltacr(i,1:3) = ( rho_b*(time_prof(i+1)-time_prof(i)) )';   
end

