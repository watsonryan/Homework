function [lat_prof,lon_prof,height_prof] = profconv(profile,orgllh)
%PROFCONV		Convert vehicle trajectory profile from local-level
%               coordinates to earth-referenced coordinates.
%       
%	[lat_prof,lon_prof,height_prof] = profconv(profile,orgllh)
%
%   INPUTS
%       profile = local-level flight profile (see PROGEN or PROGENF16)
%       orgllh = vector of x,y,z ECEF coordinates of the local-level
%                coordinate frame origin
%
%   OUTPUTS
%       lat_prof = vector of latitudes for each waypoint in the path (radians)
%       lon_prof = vector of longitudes for each waypoint in the path (radians)
%       height_prof = vector of altitudes for each waypoint in the path
%       (meters)
%

%	M. & S. Braasch 03-2005
%	Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

orgshift = [0 0 0];
llhorg = orgllh;
pos = profile(:,1:3);
L = size(pos,1);
h = waitbar(0,'ENU to LLH Trajectory Conversion');
for i = 1:L,
    xyzorg = llh2xyz(llhorg);
    xyz = enu2xyz(pos(i,1:3)-orgshift,xyzorg);
    llh = xyz2llh(xyz);
    lat_prof(i) = llh(1);
    lon_prof(i) = llh(2);
    height_prof(i) = llh(3);
    llhorg = llh;
    orgshift = pos(i,1:3);
    if rem(i,1000) == 0,
       waitbar(i/L,h)
    end
end
%  Note: cannot simply allow height to be defined as above since
%        there will be a slight height growth due to the conversion
%        from local-level-tangent plane coordinates (ENU) to ellipsoidal
%        coordinates (LLH)
height_prof = height_prof(1) + ( profile(:,3) - profile(1,3) );
close(h)
