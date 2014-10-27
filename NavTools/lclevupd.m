function [DCM_ll_I, omega_el_L, omega_ie_L] = lclevupd(lat1,lat2,vx1,vx2,vy1,vy2,...
                  height1,height2,td12,tdex,tdint,DCMel,vertmech,procflg,earthflg)
%LCLEVUPD      Compute the direction cosine matrix relating local-level-frame 
%              motion (relative to the inertial frame) over an update interval.
%
%              Note that two effects are taken into account here.  One is the
%              motion of the vehicle relative to the earth (i.e., craft rate).
%              Second is the motion of the earth (and thus the local level
%              frame also) relative to the inertial frame
%
%	[DCM_ll_I,omega_el_L,omega_ie_L] = lclevupd(lat1,lat2,vx1,vx2,vy1,vy2,...
%           height1,height2,td12,tdex,tdint,DCMel,vertmech,procflg,earthflg)
%
%   INPUTS
%       NOTE: The '1' and '2' in the input arguments refer to time indices.
%             For example, lat1 is the latitude at time index 1 and lat2
%             is the latitude at time index 2.  The convention here is
%             that time index 1 is earlier than 2.  Usually, these indices
%             refer to the previous two position/velocity updates
%       latx = geodetic latitude in radians
%       vxx = local-level-frame x-coordinate velocity in m/s
%       vyx = local-level-frame y-coordinate velocity in m/s
%       heightx = vehicle height above the reference ellipsoid in meters
%       td12 = time difference (in seconds) between time indices 1 and 2
%              (this is a positive number; i.e., td12 = time2 - time1)
%       tdex = time difference between time index 2 and the median time
%              of the nav-frame update interval (i.e., extrapolated time point)
%              (this is a positive number;
%               i.e., tdex = extrapolated_time_point - time2)
%       tdint = navigation frame update interval (in seconds)
%       DCMel = 3x3 direction cosine matrix providing the
%             transformation from the earth frame
%             to the local-level (ENU) frame
%       vertmech = argument to specify vertical mechanization
%                  0 = north-pointing (default)
%                  1 = free azimuth (vertical craft rate set to zero)
%                  2 = foucault (spatial rate set to zero)
%                  3 = unipolar (wander angle rate set to +/- longitude
%                                rate; - for northern hemisphere;
%                                + for southern hemisphere)
%        procflg = processing flag; 0=first-order solution; 1=exact solution
%        earthflg = earth shape flag
%                   0 = spherical earth; 1 = WGS-84 ellipsoid (see CRAFRATE)
%
%   OUTPUTS
%       DCM_ll_I = the direction cosine matrix relating the local-level frame (ENU)
%             at the beginning of the update interval to the local-level frame
%             at the end of the update interval (relative to the inertial frame)
%       omega_el_L = craft rate (a.k.a., transport rate) vector; this is the
%             rotation rate of the local-level-frame relative to the earth frame
%             and it is expressed in local-level-frame coordinates (rad/sec)
%

%	M. & S. Braasch 4-98
%	Copyright (c) 1997-98 by GPSoft
%	All Rights Reserved.
%

if nargin<15,error('insufficient number of input arguments'),end
if (earthflg ~= 0) && (earthflg ~= 1), error('EARTHFLG not specified correctly'),end

% lat_ex = extrapol(lat1,lat2,td12,tdex);
% vx_ex = extrapol(vx1,vx2,td12,tdex);
% vy_ex = extrapol(vy1,vy2,td12,tdex);
% height_ex = extrapol(height1,height2,td12,tdex);
lat_int = interpol(lat1,lat2,td12,tdex);
vx_int = interpol(vx1,vx2,td12,tdex);
vy_int = interpol(vy1,vy2,td12,tdex);
height_int = interpol(height1,height2,td12,tdex);

% omega_el_L = crafrate(lat_ex,vx_ex,vy_ex,height_ex,DCMel,earthflg,vertmech);
omega_el_L = crafrate(lat_int,vx_int,vy_int,height_int,DCMel,earthflg,vertmech);
omega_ie_E = [0 0 7.292115e-5]';
omega_ie_L = DCMel*omega_ie_E;
omega_il_L = omega_ie_L + omega_el_L;

ang_vect = omega_il_L*tdint;
S = skewsymm(ang_vect);
if procflg == 0,
   DCM_ll_I = eye(3) - S;   % First order approximation
else
   magn = norm(ang_vect);
   if magn == 0,
      A = eye(3);
   else                   % Exact solution
      A = eye(3) - (sin(magn)/magn)*S + ( (1-cos(magn))/magn^2 )*S*S;
   end
   DCM_ll_I = A;
end
