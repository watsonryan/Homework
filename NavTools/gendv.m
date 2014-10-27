function deltav_b = gendv(vel_prof_L,DCMnb_prof)
%GENDV		Generate that component of the Delta-V corresponding
%           to earth-referenced velocity changes.  Coriolis and
%           gravity components are not computed here 
%           (see GENDVCOR)
%       
%	deltav_b = gendv(vel_prof_L,DCMnb_prof)
%
%   INPUTS
%       vel_prof_L = velocity profile over time
%          vel_prof(i,1:3) = elements of the i-th velocity vector
%                            in local-level frame (ENU if north slave)
%
%       DCMnb_prof = profile of direction cosine elements over time
%          DCMnb_prof(i,1:9) = elements of the i-th direction cosine matrix
%                            (DCM) for vehicle attitude (navigation-to-
%                            body); 1 = DCM(1,1),
%                            2 = DCM(1,2), 3 = DCM(1,3),
%                            4 = DCM(2,1), et cetera
%
%   OUTPUTS
%       deltav_b = profile of delta-V's in body frame (Nose-Rt.Wing-Down)
%

%	M. & S. Braasch 1-98
%	Copyright (c) 1998 by GPSoft
%	All Rights Reserved.
%

if nargin<1,error('insufficient number of input arguments'),end

C = [0 1 0; 1 0 0; 0 0 -1];
deltv_L = diff(vel_prof_L);
for k = 2:size(DCMnb_prof,1),
   dcmnb=[DCMnb_prof(k,1:3); DCMnb_prof(k,4:6); DCMnb_prof(k,7:9)];
   deltav_b(k-1,1:3) = ( dcmnb*(C*deltv_L(k-1,1:3)') )';
end
