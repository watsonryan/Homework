function [quanew,r] = quaupdat(quaold,ang_vect)
%QUAUPDAT		Update quaternion body-to-nav attitude vector 
%       
%	quanew = quaupdat(quaold,ang_vect)
%     or
%  [quanew,r] = quaupdat(quaold,ang_vect)
%
%   INPUTS
%       quaold  = Body-to-nav quaternion (in nav coordinates)
%                 valid at the previous update time
%
%       ang_vect = angle vector given by the integral of the
%                  body turn-rate vector (i.e., gyro outputs)
%          ang_vect(1) = x-component (nose positive);
%          ang_vect(2) = y-component (right wing positive);
%          ang_vect(3) = z-component (down positive);
%
%   OUTPUTS
%       quanew = updated quaternion
%
%       r = quaternion representing the transformation of body
%           axes at the previous time to the body axes at the
%           current time
%

%	M. & S. Braasch 1-98
%	Copyright (c) 1997-98 by GPSoft
%	All Rights Reserved.
%

if nargin<2,error('insufficient number of input arguments'),end

magn = norm(ang_vect);
if magn == 0,
   r = [1 0 0 0];
else,
   ac = cos(magn/2);
   as = ( sin(magn/2) )/magn;
   r = [ac as*ang_vect(1) as*ang_vect(2) as*ang_vect(3)];
end
quanew = quamult(quaold,r);
