function DCMnb_prof = dcmnbgen(tc_prof)
%DCMNBGEN   Function to generate a direction cosine matrix
%           (nav-frame to body-frame) profile for a constant
%           altitude, constant velocity flight path.
%
%           Although the emulated flight path is essentially
%           straight-and-level, the DCMnb profile accounts
%           for the course changes (as given by TC_PROF) which
%           are required at each waypoint on the great
%           circle route.
%       
%	DCMnb_prof = dcmnbgen(tc_prof)
%
%   INPUTS
%       tc_prof = profile of true course from each waypoint
%                 on the great circle route (radians relative
%                 to north; east directions are positive)
%
%   OUTPUTS
%       DCMnb_prof(i,1:9) = elements of the direction cosine matrix
%                         (DCM) for vehicle attitude (nav to body); 
%                         1 = DCM(1,1), 2 = DCM(1,2), 3 = DCM(1,3),
%                         4 = DCM(2,1), et cetera
%

%	M. & S. Braasch 6-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<1,error('insufficient number of input arguments'),end

for i = 1:max(size(tc_prof)),

   dcmnb = eulr2dcm([0 0 tc_prof(i)]);

   DCMnb_prof(i,1) = dcmnb(1,1);
   DCMnb_prof(i,2) = dcmnb(1,2);
   DCMnb_prof(i,3) = dcmnb(1,3);
   DCMnb_prof(i,4) = dcmnb(2,1);
   DCMnb_prof(i,5) = dcmnb(2,2);
   DCMnb_prof(i,6) = dcmnb(2,3);
   DCMnb_prof(i,7) = dcmnb(3,1);
   DCMnb_prof(i,8) = dcmnb(3,2);
   DCMnb_prof(i,9) = dcmnb(3,3);
   
end
