function [profile,errflg] = ...
   proroll(initpos,initvel,initacc,initdcm,rollv,deltat)
%PROROLL        Flight profile sub-generator for a transition
%               between banked flight and straight and
%               level flight.  Local-level (i.e., East-
%               North-Up coordinates) version suitable for 
%               short-distance, short-duration flights.
%       
%  [profile,errflg] = ...
%          proroll(initpos,initvel,initacc,initdcm,rollv,deltat)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (3 ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (3 ENU cartesian 
%                 coordinates) (m/s)
%       initacc = initial acceleration vector (this will be zeros
%                 since it is assumed that the vehicle begins the
%                 manuever in a constant velocity condition; it is
%                 included here for consistency with other
%                 progenX functions (m/s^2)
%       initdcm = initial direction cosine matrix for
%                 vehicle attitude (navigation to body
%                 frame) (3x3 matrix)
%       rollv:  rollv(1) = target roll angle (degrees);
%               rollv(2) = roll rate (degrees/sec)
%       deltat = time increment in seconds
%
%   OUTPUTS
%       profile = flight profile
%          profile(i,1:3) = ENU path generated; 1=x, 2=y, 3=z
%          profile(i,4:6) = ENU velocity; 4 = x-velocity,
%                           5 = y-velocity, 6 = z-velocity 
%          profile(i,7:9) = ENU acceleration; 7 = x-acceleration,
%                           8 = y-acceleration, 9 = z-acceleration 
%          profile(i,10:18) = elements of the direction cosine matrix
%                            (DCM) for vehicle attitude; 10 = DCM(1,1),
%                            11 = DCM(1,2), 12 = DCM(1,3),
%                            13 = DCM(2,1), et cetera
%        errflg = 0 if there is no error condition
%               = 1 if target roll angle is nearly identical to the
%                   input roll angle in which case the function exits
%                   without adding to the flight profile
%
%   NOTES
%       - Heading is kept constant as roll is changed
%       - This program is primarily used as a transition to or from
%         straight-and-level flight when making turns (see also
%         PROROLL.M)
%

%	M. & S. Braasch 3-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<6,deltat=1;end
if nargin<5,error('insufficient number of input arguments'),end
[m,n]=size(initpos); if m>n, pos=initpos'; else, pos=initpos; end
[m,n]=size(initvel); if m>n, vel=initvel'; else, vel=initvel; end

errflg = 0;
dcmbn=initdcm';
eulvect=dcm2eulr(dcmbn);
phi=eulvect(1); theta=eulvect(2); psi=eulvect(3);

prevpos = pos;

   % Determine number of points in roll transition
   if (abs(rollv(1)-phi)<1e-6),
      npts_phi = 0;
      errflg = 1;
   else,
      phi_diff = rollv(1)*pi/180 - phi;
      if rollv(2) == 0,
         error('Cannot execute roll manuever with roll rate set to 0')
      end
      npts_phi = round((abs(phi_diff)/(rollv(2)*pi/180))/deltat);
      phi_inc = phi_diff/npts_phi;
   end,
   if npts_phi > 1e4, 
      error('input error: too many points required for roll manuever')
   end
   
   i = 0;
   for kk = 1:npts_phi,
       i=i+1;

       profile(i,4:6) = vel;
       profile(i,7:9) = [0 0 0];
       profile(i,1:3) = prevpos + vel*deltat;
       prevpos = profile(i,1:3);
       
       phi = phi + phi_inc;
       dcm = eulr2dcm([phi theta psi]);
       
       profile(i,10) = dcm(1,1);
       profile(i,11) = dcm(1,2);
       profile(i,12) = dcm(1,3);
       profile(i,13) = dcm(2,1);
       profile(i,14) = dcm(2,2);
       profile(i,15) = dcm(2,3);
       profile(i,16) = dcm(3,1);
       profile(i,17) = dcm(3,2);
       profile(i,18) = dcm(3,3);
   end,
