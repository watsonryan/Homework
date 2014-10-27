function vel_new = velupdat_tmp(vel_old2,vel_old1,td12,tdex,del_Vl,...
                             omega_el_l,DCMel,g,procflg,deltat)
%VELUPDAT		Use delta-V's, et al, to update velocity. 
%       
%	vel_new = velupdat(vel_old2,vel_old1,td12,tdex,del_Vl,...
%                      omega_el_l,DCMel,g,procflg,deltat)
%
%   INPUTS
%       vel_old2 = user velocity vector (earth referenced velocity
%                 expressed in local-level frame coordinates) valid at
%                 end of the previous update interval
%       vel_old1 = user velocity vector TD12 seconds before the previous update
%       td12 = time difference (seconds) between times of validity of
%              vel_old1 and vel_old2
%       tdex = time difference between time index 2 and the median time
%              of the velocity update interval (i.e., extrapolated time point)
%              (this is a positive number;
%               i.e., tdex = extrapolated_time_point - time2)
%       del_Vl = delta-V vector expressed in the local-level frame
%       omega_el_l = craft (i.e. transport) rate expressed in the 
%                    local-level coordinates
%       DCMel = direction cosine matrix relating the earth-frame to
%               the local-level frame
%       g = value of computed gravity for the current user position
%       procflg = processing flag; '1' causes the velocity to
%                be updated without the transport rate and earth rate
%                corrections; '0' updates velocity with all corrections
%                (thus, this parameter is normally set to '0')
%       deltat = update interval (seconds)
%
%   OUTPUTS
%       vel_new = updated user velocity vector (earth referenced velocity
%                               expressed in local-level frame coordinates)
%           vel_new(1,1) = x-coordinate velocity (m/sec)
%           vel_new(2,1) = y-coordinate velocity (m/sec)
%           vel_new(3,1) = z-coordinate velocity (m/sec)

%   REFERENCES
%       Kayton, M. and W. Fried, AVIONICS NAVIGATION SYSTEMS, 2nd edition,
%       John Wiley & Sons, New York, 1997.
%
%       Titterton, D. and J. Weston, STRAPDOWN INERTIAL NAVIGATION
%       TECHNOLOGY, Peter Peregrinus Ltd. on behalf of the Institution
%       of Electrical Engineers, London, 1997.
%
%	M. & S. Braasch 4-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<10,error('insufficient number of input arguments'),end

g_vect_L = [0 0 -g]';       
omega_ie_e = 0*[0 0 7.292115e-5]';
omega_ie_L = DCMel*omega_ie_e;

for i = 1:3,
   v_ex(i) = extrapol(vel_old1(i),vel_old2(i),td12,tdex);
end

S = skewsymm(v_ex);
if procflg == 0,
   vel_new = vel_old2' + del_Vl + ( S*(omega_el_l  + 2*omega_ie_L) + g_vect_L )*deltat;
else,
   vel_new = vel_old2' + del_Vl + ( g_vect_L )*deltat;
end
