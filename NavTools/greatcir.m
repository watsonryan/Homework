function [lat_prof,lon_prof,tc_prof,totdist] = ...
   greatcir(latlonstart,latlonend,height,earthflg,dinc)
%GREATCIR   Function to calculate a great circle path
%
%  [lat_prof,lon_prof,tc_prof,totdist] = ...
%       greatcir(latlonstart,latlonend,height,earthflg,dinc)
%
%  INPUTS
%         LATLONSTART = two element array with the latitude and longitude
%                       (in radians) of the start point
%         LATLONEND = two element array with the latitude and longitude
%                       (in radians) of the end point
%         HEIGHT = altitude (in meters) for constant altitude flight path
%         EARTHFLG = Earth shape flag
%                    0 = spherical earth with radius 6378137 meters
%                    1 = ellipsoidal earth (WGS-84 parameters)
%         DINC = distance between waypoints (nautical miles).
%                This argument is optional and is set to TOTDIST/100 by default
%
%  OUTPUTS
%         LAT_PROF = vector of latitudes for each waypoint in the path (radians)
%         LON_PROF = vector of longitudes for each waypoint in the path (radians)
%         TC_PROF = vector of true course from each waypoint to the next (radians)
%         TOTDIST = total great circle path distance from origin to destination
%                   (nautical miles)
%
%  NOTE:  NEITHER THE START POINT (ORIGIN) NOR THE END POINT (DESTINATION)
%         CAN BE AT THE EARTH'S POLES
%
%         Also, paths which go over a pole exactly will have a slight discontinuity
%
%  REFERENCES
%         Aviation Formulary V1.21 by Ed Williams, 
%                  http://www.best.com/~williams/avform.htm
%
%         Siouris, G., AEROSPACE AVIONICS SYSTEMS - A MODERN SYNTHESIS,
%         Academic Press, Inc., San Diego, 1993.

%  
%	M. & S. Braasch 6-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%
if nargin<5, distflg = 0; else distflg = 1; end
if nargin<4,error('insufficient number of input arguments'),end
if (earthflg ~= 0) && (earthflg ~= 1), error('EARTHFLG not specified correctly'),end

lat_start = latlonstart(1); lon_start = latlonstart(2);
lat_end = latlonend(1); lon_end = latlonend(2);
lat1 = lat_start; lon1 = lon_start;
lat2 = lat_end; lon2 = lon_end;
%
if earthflg == 0,
   Ravg_m = 6378137 + height;
else
   ae = 6378137;
   e = 0.0818191908426;
   x1 = ( ae/(sqrt(1-e*e*sin(lat1)*sin(lat1))) + height )*cos(lat1);
   y1 = ( ae*(1-e*e)/(sqrt(1-e*e*sin(lat1)*sin(lat1))) + height )*sin(lat1);
   x2 = ( ae/(sqrt(1-e*e*sin(lat2)*sin(lat2))) + height )*cos(lat2);
   y2 = ( ae*(1-e*e)/(sqrt(1-e*e*sin(lat2)*sin(lat2))) + height )*sin(lat2);
   R1 = norm([x1 y1]);
   R2 = norm([x2 y2]);
   Ravg_m = mean([R1 R2]);
end
Ravg_nmi = (Ravg_m/0.3048)/6076;
%
d_rad = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2));
d_nm = d_rad*Ravg_nmi;
%
if distflg == 0,
   d_inc_nm = d_nm/100;
else
   d_inc_nm = dinc;
end
d_inc_rad = d_inc_nm/Ravg_nmi;

N = ceil(d_nm/d_inc_nm);
for i = 1:N,
   d_rad = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2));
   tc = acos((sin(lat2)-sin(lat1)*cos(d_rad))/(sin(d_rad)*cos(lat1)));
   if sin(lon2-lon1) < 0,
      tc = 2*pi-tc;
   end
   if ~isreal(tc), tc = real(tc); end
   if tc > (2*pi - 1e-4), tc = 0; end
   tc_prof(i) = tc;

   lat(i) = asin(sin(lat1)*cos(d_inc_rad)+cos(lat1)*sin(d_inc_rad)*cos(tc));
   if abs(cos(lat(i))) < 3*eps,
      lon(i) = lon1;
   else
      atmp = sin(tc)*sin(d_inc_rad)/cos(lat(i));
      if atmp < -1,
         lon(i) = mod(lon1+asin(-1)+pi,2*pi)-pi;
      elseif atmp > 1,
         lon(i) = mod(lon1+asin(1)+pi,2*pi)-pi;
      else
         lon(i) = mod(lon1+asin(sin(tc)*sin(d_inc_rad)/cos(lat(i)))+pi,2*pi)-pi;
      end
   end
   lat1 = lat(i);
   lon1 = lon(i);
end
lat_prof = lat;
lon_prof = lon;
totdist = d_nm;
