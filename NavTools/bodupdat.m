function [DCMbn_new,DCMbb] = bodupdat(DCMbn_old,ang_vect)
%BODUPDAT		Update the direction cosine matrix
%              for body motion (relative to inertial space).
%              The function is thus acting upon the strapdown
%              gyro outputs.
%       
%	DCMbn_new = bodupdat(DCMbn_old,ang_vect)
%     or
%  [DCMbn_new,DCMbb] = bodupdat(DCMbn_old,ang_vect)
%
%   INPUTS
%       DCMbn_old  = current direction cosine matrix providing the
%                    transformation from body to nav coordinates
%
%       ang_vect = incremental integral of body angular
%                  rate vector; in the absence of coning
%                  (i.e., angular rate vector is constant
%                  over the integration interval), this
%                  is the output of the rate-integrating
%                  gyros.
%          ang_vect(1) = x-component (roll);
%          ang_vect(2) = y-component (pitch);
%          ang_vect(3) = z-component (yaw);
%
%   OUTPUTS
%       DCMbn_new = updated direction cosine matrix relating the
%                   current body-frame to the previous
%                   navigation-frame.  Note that the update of
%                   the navigation frame (which is related to
%                   the local-level frame) is accomplished in
%                   function LCLEVUPD
%
%       DCMbb = direction cosine matrix providing the
%               transformation from the body coordinates at
%               time k+1 to the body coordinates at time k
%

%	M. & S. Braasch 4-98
%	Copyright (c) 1997-98 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<2,error('insufficient number of input arguments'),end

S = skewsymm(ang_vect);
magn = norm(ang_vect);
if magn == 0,
   DCMbb = eye(3);
else
   DCMbb = eye(3) + (sin(magn)/magn)*S + ( (1-cos(magn))/magn^2 )*S*S;
end
DCMbn_new = DCMbn_old*DCMbb;
