function F_ins = F_ins_gen(lat,del_V_L,td,omega_el_L,omega_ie_L,grav,height)
%F_INS_GEN		Generate the 9x9 continuous-time system dynamics
%               matrix for the inertial error model with PSI-angle
%               attitude error representation.
%
%               Position and velocity errors are expressed in the local
%               level frame known as the L-frame with x=east; y=north and
%               z=up
%       
%	F_ins = F_INS_GEN(LAT,DEL_V_L,TD,OMEGA_EL_L,OMEGA_IE_L,GRAV,HEIGHT)
%
%   INPUTS
%       LAT = latitude (in radians)
%       DEL_V_L = delta-V vector converted to the L-frame  (m/s)
%       TD = delta-V integration interval in seconds
%       OMEGA_EL_L = transport (or craft) rate vector in the L-frame
%                                                           (rad/s)
%       OMEGA_IE_L = earth-rate vector in the L-frame (rad/s)
%       GRAV = local value of gravity (in meters/sec)
%       HEIGHT = height above the WGS-84 ellipsoid
%
%   OUTPUTS
%       F_ins = 9x9 continuous-time system dynamics matrix for
%               the inertial error model
%

%  REFERENCES
%
%   Huddle, J. R., "Inertial Navigation System Error Model Considerations
%   in Kalman Filtering Applications," CONTROL AND DYNAMIC SYSTEMS, Vol.
%   20, Academic Press, 1983.
%
%   Pue, A., "Integration of GPS with Inertial Navigation Systems,"
%   short course notes, Navtech Seminars, Springfield, VA, 2003.
%
%	M. Braasch    September 2005
%	Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<7,error('insufficient number of input arguments'),end

[rm,rp] = radicurv(lat);  radius_e = sqrt(rm*rp);
accel_vect_L = del_V_L*(1/td);
F11 = -1*skewsymm(omega_el_L);    % L-frame
F12 = eye(3);
F13 = zeros(3,3);
F21=eye(3); F21(1,1)=-grav/radius_e; 
F21(2,2)=-grav/radius_e; F21(3,3)=2*grav/(radius_e+height);
F22 = -1*skewsymm(2*omega_ie_L + omega_el_L);   % L-frame
F23 = skewsymm(accel_vect_L);
F31 = zeros(3,3);
F32 = zeros(3,3);
F33 = (-1)*skewsymm(omega_ie_L + omega_el_L);   % L-frame
F_ins = [F11 F12 F13; F21 F22 F23; F31 F32 F33];

