function rho = crafrate(lat,vx,vy,height,DCMel,earthflg,vertmech)
%CRAFRATE		Calculate craft (a.k.a. transport) rate. 
%       
%	rho = crafrate(lat,vx,vy,height,DCMel,earthflg,vertmech)
%
%   INPUTS
%       lat = geodetic latitude in radians
%       vx = local-level-frame x-coordinate velocity in m/s
%       vy = local-level-frame y-coordinate velocity in m/s
%       height = vehicle height above the reference ellipsoid in meters
%       DCMel = direction cosine matrix relating the earth-frame to
%               the local-level frame
%       earthflg = earth shape flag
%                  0 = spherical earth with radius 6378137 meters
%                  1 = ellipsoidal earth (WGS-84)
%       vertmech = Optional argument specifying vertical mechanization
%                  0 = north-pointing (default)
%                  1 = free azimuth (vertical craft rate set to zero)
%                  2 = foucault (vertical component of spatial rate set to zero)
%                  3 = unipolar (wander angle rate set to +/- longitude
%                                rate; - for northern hemisphere;
%                                + for southern hemisphere)
%
%   OUTPUTS
%       rho = craft rate vector
%           rho(1,1) = local-level-frame x-coordinate angular rate (radians/sec)
%           rho(2,1) = local-level-frame y-coordinate angular rate (radians/sec)
%           rho(3,1) = local-level-frame z-coordinate angular rate (radians/sec)
%
%   NOTES
%       For the north-pointing case, rho(1) is east craft rate and
%       rho(2) is north craft rate.  For all cases, rho(3) is the
%       (positive) up craft rate.
%
%       For the wander azimuth cases, first-order approximations are
%       implemented as per Brockstein and Kouba.
%
%   REFERENCES
%       Kayton, M. and W. Fried, AVIONICS NAVIGATION SYSTEMS, 2nd edition,
%       John Wiley & Sons, New York, 1997.
%
%       Titterton, D. and J. Weston, STRAPDOWN INERTIAL NAVIGATION
%       TECHNOLOGY, Peter Peregrinus Ltd. on behalf of the Institution
%       of Electrical Engineers, London, 1997.
%
%       Brockstein, A. and J. Kouba, "Derivation of Free Inertial, General
%       Wander Azimuth Mechanization Equations," Litton Systems, Inc., 
%       Guidance and Control Systems Division, Woodland Hills, California,
%       June 1969, Revised June 1981.

%	M. & S. Braasch 8-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<7, vertmech = 0; end
if nargin<6,error('insufficient number of input arguments'),end
if (earthflg ~= 0) && (earthflg ~= 1), error('EARTHFLG not specified correctly'),end

if earthflg == 0,
   rm = 6378137;
   rp = rm;
else
   [rm,rp] = radicurv(lat);
end

if vertmech == 0,
   Ve = vx;
   Vn = vy;
   rho(1,1) = -Vn/(rm + height);
   rho(2,1) = Ve/(rp + height);
   rho(3,1) = Ve*tan(lat)/(rp + height);
else
   h = height;
   f = 1/298.257223563;
   [b,a] = radicurv(0);
   DCMelYx = DCMel(1,3);
   DCMelYy = DCMel(2,3);
   DCMelYz = DCMel(3,3);

   rho(1,1) = -(vy/a)*(1 - h/a - f*(1-3*DCMelYy*DCMelYy-DCMelYx*DCMelYx)) ...
      - (vx/a)*(2*f*DCMelYx*DCMelYy);
%
   rho(2,1) = (vx/a)*(1 - h/a - f*(1-3*DCMelYx*DCMelYx-DCMelYy*DCMelYy)) ...
      + (vy/a)*(2*f*DCMelYx*DCMelYy);

   if vertmech == 1,
     rho(3,1) = 0;
   end

   if vertmech == 2,
      OMEGA = 7.292115e-5;   %Earth's rate in rad/s
      rho(3,1) = -OMEGA*DCMelYz;
   end

   if vertmech == 3,
      num = -( rho(2,1)*DCMelYy + rho(1,1)*DCMelYx );
      if lat >= 0,
         den = DCMelYz + 1;
      else
         den = DCMelYz - 1;
      end
      rho(3,1) = num/den;
   end
   
end
