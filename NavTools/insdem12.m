function insdem12(latstart,lonstart,latend,lonend)
%INSDEM12      Function to illustrate the use of the routines
%              which generate and display great circle paths
%
%    INSDEM12(latstart,lonstart,latend,lonend)
%
%  INPUTS
%      LATSTART: latitude of starting point in degrees
%      LONSTART: longitude of starting point in degrees
%      LATEND: latitude of ending point in degrees
%      LONEND: longitude of ending point in degrees
%
%      The input arguments are optional.  The default
%      starting and ending points are New York and
%      Istanbul
%
%      By convention, north latitudes and east longitudes
%      are positive
%

%	M. & S. Braasch 7-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

%
if nargin < 4,
   flag = 1;
else
   flag = 0;
end
%
d2r = pi/180;
%
JFK_deg = [40+38/60 -(73+47/60)];   % New York (JFK Airport)
JFK_rad = JFK_deg*d2r;
IST_deg = [40.983 28.817];         % Istanbul
IST_rad = IST_deg*d2r;

if flag == 1,
   latlonstart = [JFK_rad(1) JFK_rad(2)];
   latlonend = [IST_rad(1) IST_rad(2)];
else
   latlonstart = [latstart*d2r lonstart*d2r];
   latlonend = [latend*d2r lonend*d2r];
end
%
[lat_vec, lon_vec, tc_vec] = greatcir(latlonstart,latlonend,0,1);
pathplot(lat_vec,lon_vec)
