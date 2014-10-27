function deltaer = earthrot(timevec,DCMel_prof,DCMnb_prof)
%EARTHROT		Generate ideal earth rotation measurements as a function of
%              user latitude, user longitude, measurement time interval
%              and user-body attitude
%       
%	deltaer = earthrot(timevec,DCMel_prof,DCMnb_prof)
%
%   INPUTS
%       timevec = vector of time points associated with the user's
%                 latitude and longitude as specified in LAT_PROF 
%                 and LONG_PROF (seconds)
%       DCMel_prof = profile of earth-frame to local-level-frame (ENU) direction 
%                    cosine elements over time
%          DCMel_prof(i,1:9) = elements of the i-th direction cosine matrix
%                            (DCM); 1 = DCM(1,1),
%                            2 = DCM(1,2), 3 = DCM(1,3),
%                            4 = DCM(2,1), et cetera
%       DCMnb_prof = profile of nav-frame (NED) to body-frame direction 
%                    cosine elements over time
%          DCMnb_prof(i,1:9) = elements of the i-th direction cosine matrix
%                            (DCM) for vehicle attitude; 1 = DCM(1,1),
%                            2 = DCM(1,2), 3 = DCM(1,3),
%                            4 = DCM(2,1), et cetera
%
%   OUTPUTS
%       deltaer = vector of ideal earth rotation measurements (radians);
%                 If the input vectors are N points long, DELTAER is
%                 then N-1 points long.  DELTAER(1,:) corresponds to the
%                 time specified by TIMEVEC(2).  DELTAER(2,:) corresponds
%                 to TIMEVEC(3) and so on.  For a given index 'i', 
%                 DELTAER(i,1:3) are the three components of the earth
%                 rotation vector expressed in the body-frame.
%

%	M. & S. Braasch 5-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<3,error('insufficient number of input arguments'),end

omega_ie_E = [0 0 7.292115e-5]';
C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

deltatvec = diff(timevec);

for k = 2:max(size(timevec)),
  DCMnb2=[DCMnb_prof(k,1:3); DCMnb_prof(k,4:6); DCMnb_prof(k,7:9)];
  DCMnb1=[DCMnb_prof(k-1,1:3); DCMnb_prof(k-1,4:6); DCMnb_prof(k-1,7:9)];
  DCMnbavg = 0.5*( DCMnb2 + DCMnb1 );
  DCMel2=[DCMel_prof(k,1:3); DCMel_prof(k,4:6); DCMel_prof(k,7:9)];
  DCMel1=[DCMel_prof(k-1,1:3); DCMel_prof(k-1,4:6); DCMel_prof(k-1,7:9)];
  DCMelavg = 0.5*( DCMel2 + DCMel1 );
  tmp1 = omega_ie_E*deltatvec(k-1);     % integral of earth rotation over time interval
  tmp2 = DCMelavg*tmp1;             % convert from earth-frame to local-level-frame (ENU)
  tmp3 = C*tmp2;                 % convert from local-level-fram (ENU) 
  %                              % to the navigation frame (NED)
  deltaer(k-1,1:3) = (DCMnbavg*tmp3)'; % convert from nav-frame to body-frame
end
