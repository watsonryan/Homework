function pathplot(lat_vec,lon_vec)
%PATHPLOT   Function to plot a sequence of waypoints on the globe
%
%  pathplot(lat_vec,lon_vec)
%
%  INPUTS
%         LAT_VEC = array of latitudes (in radians) of the waypoints
%         LON_VEC = array of longitudes (in radians) of the waypoints

%  
%	M. & S. Braasch 6-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

hold on
plotwrld
N = max(size(lat_vec));
lon = lon_vec;
lat = lat_vec;
lat_start = lat(1);  lat_end = lat(N);
lon_start = lon(1);  lon_end = lon(N);
for j = 1:N,
   [xp(j),yp(j),zp(j)] = sph2cart(lon(j),lat(j),1);
end

h=plot3(xp,yp,zp,'r.', ...
	'MarkerSize',10, ...
   'EraseMode','none');
%drawnow

K = 89.99*pi/180;
if (abs(lon_start)>K)&&(abs(lon_end)>K)
   if lon_start < 0,
      lona=lon_start+2*pi; 
   else
      lona=lon_start; 
   end
   if lon_end < 0, 
      lonb=lon_end+2*pi; 
   else
      lonb=lon_end; 
   end
   flag = 0;
elseif (lon_start<0)&&(lon_end>K)
   lona = lon_start;
   lonb = lon_end;
   flag = 1;
elseif (lon_start>K)&&(lon_end<0)
   lona = lon_start;
   lonb = lon_end;
   flag = 1;
else
   lona = lon_start;
   lonb = lon_end;
   flag = 0;
end


azview = (180/pi)*(pi/2 + mean([lona lonb]));
if flag == 1, azview = azview - 180; end
elview = (180/pi)*mean([lat_start lat_end]);
view(azview,elview)
drawnow
hold off
