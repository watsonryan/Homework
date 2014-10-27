function DCMel = llw2dcm(llw_vect)
%LLW2DCM       Convert from latitude-longitude-wander_angle
%              to the direction cosine matrix relating the
%              earth frame to the local-level (ENU) frame.
%       
%	DCMel = llw2dcm(llw_vect)
%
%   INPUTS
%       llw_vect(1) = latitude in radians 
%
%       llw_vect(2) = longitude in radians 
%
%       llw_vect(3) = wander angle in radians 
%
%   OUTPUTS
%       DCMel = 3x3 direction cosine matrix providing the
%             transformation from the earth frame
%             to the local-level (ENU) frame
%
%              The earth frame is an earth-centered, earth-fixed
%              cartesian coordinate system with the x-axis
%              pointing out of the intersection of the prime
%              meridian and the equator; the z-axis points
%              along the north pole and the y-axis completes
%              right-hand coordinate system.
%
%              The locally level frame is a
%              cartesian coordinate system with x-y-z axes
%              pointing along East-North-Up (ENU).
%

%   REFERENCE:  Titterton, D. and J. Weston, STRAPDOWN
%               INERTIAL NAVIGATION TECHNOLOGY, Peter
%               Peregrinus Ltd. on behalf of the Institution
%               of Electrical Engineers, London, 1997.
%
%               Kayton & Fried, 1997.
%
%	M. & S. Braasch 4-98
%	Copyright (c) 1998 by GPSoft
%	All Rights Reserved.
%

  if nargin<1,error('insufficient number of input arguments'),end

  lat = -llw_vect(1); long = llw_vect(2); alpha = llw_vect(3);
  %% note the minus sign on the latitude takes care of the negative rotation

  clat = cos(lat); slat = sin(lat);
  clong = cos(long); slong = sin(long);
  calpha = cos(alpha); salpha = sin(alpha);

  C1 = [clong  slong 0; ...
        -slong clong 0; ...
         0     0   1];            % Rotation of LONG radians about z
  C2 = [clat  0  -slat; ...
          0   1     0 ; ...
        slat  0   clat];          % Rotation of LAT radians about y
  C3 = [1   0    0;   ...
        0  calpha salpha; ...
        0 -salpha calpha];        % Rotation of ALPHA radians about x

  C = C3 * C2 * C1;  % We've converted from lat-long-wander to Up-East-North
  
  Cconv = ...
     [0 1 0;
     0 0 1;
     1 0 0];   % Use this to convert from Up-East-North to East-North-Up
  
  DCMel = Cconv*C;
  