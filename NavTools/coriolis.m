function vel_cor = coriolis(vel_extrp_L,omega_en_n,DCMbn,DCMel,outflg,deltat)
%CORIOLIS		Compute the Coriolis correction relating earth-referenced
%              velocities to inertial velocities. 
%       
%	vel_cor = coriolis(vel_extr_L,omega_en_n,DCMbn,DCMel,outflg,deltat)
%
%   INPUTS
%       vel_extrp_L = user velocity vector (earth referenced velocity
%                 expressed in local-level frame coordinates) extrapolated
%                 such that it is valid at the middle of the current
%                 update interval
%       omega_en_n = craft rate vector expressed in nav coordinates
%       DCMbn = direction cosine matrix providing the transformation
%               from body to nav coordinates
%       DCMel = direction cosine matrix providing the transformation
%               from earth to local-level coordinates
%       outflg = output data flag (0 = output in nav-frame coordinates
%                                  1 = output in body-frame coordinates)
%       deltat = update interval (seconds)
%
%   OUTPUTS
%       vel_cor = user velocity vector correction (earth referenced 
%                 velocity expressed in body or nav frame coordinates)
%           vel_cor(1,1) = body- or nav-frame x-coordinate velocity (m/sec)
%           vel_cor(2,1) = body- or nav-frame y-coordinate velocity (m/sec)
%           vel_cor(3,1) = body- or nav-frame z-coordinate velocity (m/sec)
%
%       The convention is as follows: vel_cor = v X (2omega + rho)
%           where v is the user velocity
%                 X denotes the cross-product operator
%                 omega is the earth-rate vector
%                 rho is the transport (a.k.a. craft) rate vector

%   REFERENCES
%       Kayton, M. and W. Fried, AVIONICS NAVIGATION SYSTEMS, 2nd edition,
%       John Wiley & Sons, New York, 1997.
%
%       Titterton, D. and J. Weston, STRAPDOWN INERTIAL NAVIGATION
%       TECHNOLOGY, Peter Peregrinus Ltd. on behalf of the Institution
%       of Electrical Engineers, London, 1997.
%
%	M. & S. Braasch 5-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<6,error('insufficient number of input arguments'),end
C = [0 1 0; 1 0 0; 0 0 -1];      % Conversion between ENU and NED
omega_ie_e = [0 0 7.292115e-5]';
omega_ie_n = C*(DCMel*omega_ie_e);
S = skewsymm(C*vel_extrp_L');
vel_cor_n = S*(omega_en_n + 2*omega_ie_n)*deltat;
if outflg == 1,
   vel_cor = ( DCMbn' )*vel_cor_n;   % conversion to body coordinates
else
   vel_cor = vel_cor_n;
end
