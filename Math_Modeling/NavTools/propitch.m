function [profile,errflg] = ...
   propitch(initpos,initvel,initacc,initdcm,pitchv,deltat)
%PROPITCH        Flight profile sub-generator for a transition
%               between pitched flight and straight and
%               level flight.  Local-level (i.e., East-
%               North-Up coordinates) version suitable for 
%               short-distance, short-duration flights.
%       
%  [profile,errflg] = ...
%          propitch(initpos,initvel,initacc,initdcm,pitchv,deltat)
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
%       pitchv:  pitchv(1) = target pitch angle (degrees);
%                pitchv(2) = pitch rate (degrees/sec)
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
%               = 1 if pitch manuever cannot be simulated and the user should
%                 specify a smaller time step or smaller pitch rate
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

   % Determine number of points in pitch transition
   if (pitchv(1) == -999) || (abs(pitchv(1)-theta)<1e-6),
      npts_the = 0;
   else
      theta_diff = pitchv(1)*pi/180 - theta;
      if pitchv(2) == 0,
         error('Cannot execute pitch manuever with pitch rate set to 0')
      end
      r = norm(vel)/(pitchv(2)*pi/180);
      cntrpacc = ( norm(vel) )^2 /r;
      npts_the=round((abs(theta_diff)/(pitchv(2)*pi/180))/deltat);
      if npts_the<2, 
         errflg = 1;
         fprintf(1,'pitch maneuver cannot be simulated !\n')
         error('specify smaller time step or slower pitch rate')
      end
      if npts_the > 1e4,
       error('input error: too many points required for pitch manuever')
      end
      theta_inc = theta_diff/npts_the;
      velang=atan2(vel(3),norm(vel(1:2)));
      beta=velang+(pi/2)*sign(theta_diff);
      turnorg(3)=pos(1,3)+r*sin(beta);
      tmp=r*cos(beta); gamma=atan2(vel(2),vel(1));
      turnorg(1)=pos(1,1)+tmp*cos(gamma);
      turnorg(2)=pos(1,2)+tmp*sin(gamma);
      cumtheang=0;
   end,

   i = 0;
   for kk = 1:npts_the,
       i=i+1;

          cumtheang = cumtheang + theta_inc;
          profile(i,3)=turnorg(3)+r*sin(beta-pi+cumtheang);
          tmp=r*cos(beta-pi+cumtheang);
          profile(i,1)=turnorg(1)+tmp*cos(gamma);
          profile(i,2)=turnorg(2)+tmp*sin(gamma);
      
          rotmat1 = [cos(gamma) sin(gamma) 0;
                    -sin(gamma) cos(gamma) 0;
                         0          0      1];
          rotmat2 = [cos(theta_inc) 0 -sin(theta_inc);
                          0          1            0    ;
                     sin(theta_inc) 0  cos(theta_inc)];
          rotmat3 = [cos(-gamma) sin(-gamma) 0;
                    -sin(-gamma) cos(-gamma) 0;
                        0              0     1];
          vel=(rotmat3*rotmat2*rotmat1*vel')';
          profile(i,4:6) = vel;

          accdvect= turnorg(1:3) - profile(i,1:3);
          naccdvec=accdvect/norm(accdvect);
          profile(i,7:9) = cntrpacc*naccdvec;
          theta = theta + theta_inc;
 
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
