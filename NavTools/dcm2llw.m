function llw_vect = dcm2llw(DCMel)
%DCM2LLW       Conversion from direction cosine matrix
%              (earth-frame to local-level-frame) to 
%              latitude, longitude and wander angle
%       
%	 llw_vect = dcm2llw(DCMel)
%
%   INPUTS
%       DCMel = 3x3 direction cosine matrix providing the
%             transformation from the earth frame
%             to the local-level frame
%
%   OUTPUTS
%       llw_vect(1) = latitude in radians 
%
%       llw_vect(2) = longitude in radians 
%
%       llw_vect(3) = wander angle in radians 
%

%   REFERENCE:  Kayton, M. and W. Fried, editors,
%               AVIONICS NAVIGATION SYSTEMS, 2nd edition,
%               Wiley-Interscience, John Wiley & Sons,
%               New York, 1997.
%
%	M. & S. Braasch 4-98
%	Copyright (c) 1998 by GPSoft
%	All Rights Reserved.
%

  if nargin<1,error('insufficient number of input arguments'),end

  lat = asin(DCMel(3,3));
  long = atan2(DCMel(3,2),DCMel(3,1));
  alpha = atan2(DCMel(1,3),DCMel(2,3));
  
  llw_vect(1) = lat;
  llw_vect(2) = long;
  llw_vect(3) = alpha;