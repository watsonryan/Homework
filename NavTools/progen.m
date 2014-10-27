function profile = progen(initpos,initvel,initdcm,segparam)
%PROGEN		Flight profile generator.  Local-level version
%               suitable for short-distance, short-duration flights.
%               Profile can consist of:
%		- constant-velocity straight segments (level, 
%                 climbing or descending)
%		- constant-acceleration straight segments (level,
%                 climbing or descending)
%		- constant-altitude, constant-radius turns
%     - transitions between flight segments (e.g.,
%       between straight-and-level and climbing flight)
%       
%	profile = progen(initpos,initvel,initdcm,segparam)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (ENU
%                 cartesian coordinates) (m/s)
%       initdcm = initial direction cosine matrix (nav-to-body)
%                 for vehicle attitude (3x3 matrix)
%       segparam = segment and turn parameters;  N x ? matrix
%                  where ...
%                  segparam(i,1) = segment type identifier
%                       1 = straight constant-velocity
%                       2 = straight constant-acceleration
%                       3 = constant-altitude, constant-radius turn
%                       4 = transition between segments
%                  segparam(i,2) = duration (in seconds) of i-th segment;
%                                  if the segment type number is 1, the
%                                  velocity will be held constant at the
%                                  value achieved at the end of the previous
%                                  segment (except in the case of segment 
%                                  number 1 in which case the input vector
%                                  INITVEL is used);
%                                  if this segment is a turn or transition,
%                                  this parameter is ignored (i.e., treated
%                                  as a dummy argument)
%                  segparam(i,3) = total acceleration in g's for i-th 
%                                  segment; this parameter is ignored
%                                  for all segment type numbers other than 2
%                  segparam(i,4) = amount of turn (degrees); note the
%                                  direction of the turn is given by the
%                                  bank (roll) angle specified in the
%                                  direction cosine matrix (dcm); the
%                                  aircraft thus must be banked prior to
%                                  the execution of a turn; use the
%                                  transition manuever (segment type
%                                  number 4 to roll into the proper
%                                  bank angle;
%                                  for segment type numbers other than 3, 
%                                  this parameter is a dummy argument
%                  segparam(i,5) = used for transitions only; this parameter
%                                  is the desired roll (bank) angle to be 
%                                  achieved by the end of the transition
%                                  (degrees); if the roll angle is to be
%                                  left unchanged (i.e., during a pitch
%                                  only transition), this parameter should
%                                  be set to -999; for segment type numbers
%                                  other than 4, this parameter is a dummy
%                                  argument
%                  segparam(i,6) = used for transitions only; this parameter
%                                  is the roll-rate (degrees/sec) to be
%                                  used during a roll manuever; for segment
%                                  numbers other than 4, this parameter is
%                                  a dummy argument
%                  segparam(i,7) = used for transitions only; this parameter
%                                  is the desired pitch angle to be 
%                                  achieved by the end of the transition
%                                  (degrees); if the pitch angle is to be
%                                  left unchanged (i.e., during a roll
%                                  only transition), this parameter should
%                                  be set to -999; for segment type numbers
%                                  other than 4, this parameter is a dummy
%                                  argument
%                  segparam(i,8) = used for transitions only; this parameter
%                                  is the pitch-rate (degrees/sec) to be
%                                  used during a pitch manuever; for segment
%                                  numbers other than 4, this parameter is
%                                  a dummy argument
%                  segparam(i,9) = time step for i-th segment (seconds)
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
%          profile(i,19) = simulation run time (seconds)
%
%   NOTES
%       For each segment specified, all eight parameters in SEGPARAM must
%       be specified.  Dummy arguments must have a value assigned.
%
%       Roll and pitch transitions must be accomplished separately.  As a
%       result, separate segments (and, thus, rows in the SEGPARAM matrix)
%       must be specified for roll and pitch maneuvers.

%	M. & S. Braasch 12-97, Updated 8-98
%	Copyright (c) 1997-98 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<5,deltat=1;end
if nargin<4,error('insufficient number of input arguments'),end

[m,n]=size(initvel); if m>n, initvel=initvel'; end
[m,n]=size(initpos); if m>n, initpos=initpos'; end
[m,n]=size(segparam); nsegs = m;

profile(1,1:3) = initpos;
profile(1,4:6) = initvel;
profile(1,7:9) = [0 0 0];
profile(1,10:12) = initdcm(1,1:3);
profile(1,13:15) = initdcm(2,1:3);
profile(1,16:18) = initdcm(3,1:3);

for i = 1:nsegs,

   fprintf(1,' Creating profile for segment number %i \n',i)
   deltat = segparam(i,9);
   
   if i == 1,
      pos = initpos;
      vel = initvel;
      dcm = initdcm;
   else
      k = size(profile,1);
      pos = profile(k,1:3);
      vel = profile(k,4:6);
      dcm(1,1:3) = profile(k,10:12);
      dcm(2,1:3) = profile(k,13:15);
      dcm(3,1:3) = profile(k,16:18);
   end,

   if segparam(i,1) == 1,
      acc = [0 0 0];
      duration = segparam(i,2);
      profilez = prostrai(pos,vel,acc,dcm,duration,deltat);
   elseif segparam(i,1) == 2,
      eul_vect = dcm2eulr(dcm);
      z = sin(eul_vect(2));
      x = sin(eul_vect(3));
      y = cos(eul_vect(3));
      acc = segparam(i,3)*9.81*[x y z];
      duration = segparam(i,2);      
      profilez = prostrai(pos,vel,acc,dcm,duration,deltat);
   elseif segparam(i,1) == 3,
      acc = [0 0 0];
      turnamt = segparam(i,4);
      [profilez,errflg] = ...
                proturn(pos,vel,acc,dcm,turnamt,deltat);
   elseif segparam(i,1) == 4,
      acc = [0 0 0];
      if segparam(i,5) == -999,
         pitchv = [segparam(i,7) segparam(i,8)];
         [profilez,errflg] = propitch(pos,vel,acc,dcm,pitchv,deltat);
      else
         rollv = [segparam(i,5) segparam(i,6)];
         [profilez,errflg] = ...
            proroll(pos,vel,acc,dcm,rollv,deltat);
      end,
   end
   
   if i == 1,
      mm = size(profilez(:,18),1);
      profilez(:,19) = (0:deltat:deltat*(mm-1))';
      profile = profilez;
   else
      mm = size(profilez(:,18),1);
      nn = size(profile,1);
      oldtime = profile(nn,19);
      profilez(1:mm,19) = oldtime + deltat*(1:mm)';
      profile = [profile; profilez];
   end
   
end

