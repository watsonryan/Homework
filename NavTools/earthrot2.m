function deltaer = earthrot2(deltat,DCMel,DCMbn)
%EARTHROT2		Generate a single ideal earth rotation measurement as a function of
%              user latitude, user longitude, measurement time interval
%              and user-body attitude
%       
%	deltaer = earthrot2(deltat,DCMel,DCMbn)
%
%   INPUTS
%       deltat = time interval between current and previous
%                 sample points (seconds)
%       DCMel = earth-frame to local-level-frame (ENU) direction 
%                    cosine matrix
%       DCMbn = body-to-nav direction cosine matrix
%
%   OUTPUTS
%       deltaer = ideal earth rotation measurement (radians);
%                 DELTAER(1:3) are the three components of the earth
%                 rotation vector expressed in the body-frame.
%

%	M. & S. Braasch 5-98, 11-2008
%	Copyright (c) 1998, 2008 by GPSoft LLC
%	All Rights Reserved.
%
%   REVISION HISTORY
%
%   November 2008:  function created to provide a single point instead of
%   an entire profile

if nargin<3,error('insufficient number of input arguments'),end

omega_ie_E = [0 0 7.292115e-5]';
C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU
  tmp1 = omega_ie_E*deltat;     % integral of earth rotation over time interval
  tmp2 = DCMel*tmp1;             % convert from earth-frame to local-level-frame (ENU)
  tmp3 = C*tmp2;                 % convert from local-level-fram (ENU) 
  %                              % to the navigation frame (NED)
  deltaer(1,1:3) = ( (DCMbn)' * tmp3 )'; % convert from nav-frame to body-frame
end
