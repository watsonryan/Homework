function profile = progenf16(initpos,initvel,initdcm,segparam)
%PROGENF16		F-16 flight profile generator.  Local-level version
%               suitable for short-distance, short-duration flights.
%               Profile can consist of:
%		- constant-velocity straight-and-level segments
%		- constant-acceleration straight-and-level segments
%       - constant rate-of-climb segments
%		- constant-altitude, coordinated turns
%       - transitions between straight-and-level and climbing flight
%       - transitions between straight-and-level and turning flight
%       
%	profile = progenf16(initpos,initvel,initdcm,segparam)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (ENU
%                 cartesian coordinates) (m/s)
%       initdcm = initial direction cosine matrix (nav-to-body)
%                 for vehicle attitude (3x3 matrix)
%       segparam = segment and turn parameters;  N x 8 matrix
%                  where ...
%                  segparam(i,1) = segment-type identifier
%                       5 = level, straight constant-acceleration
%                           including constant-velocity as a special case)
%                       6 = transition between level and climb
%                       7 = constant rate-of-climb
%                       8 = transition between level-straight and
%                       level-turn
%                       9 = level coordinated turn
%                  segparam(i,2) = duration (in seconds) of i-th segment;
%                                  if the segment type number is 1, the
%                                  velocity will be held constant at the
%                                  value achieved at the end of the previous
%                                  segment (except in the case where this is
%                                  the first segment of the simulation in 
%                                  which case the input vector INITVEL is used);
%                                  if this segment is a turn or transition,
%                                  this parameter is ignored (i.e., treated
%                                  as a dummy argument)
%                  segparam(i,3) = final velocity in meters-per-second; this
%                                  parameter is ignored for all segment
%                                  type numbers other than 5
%                  segparam(i,4) = flight path angle (i.e., climb angle) in
%                                  degrees; this parameter is ignored for
%                                  all segment type numbers other than
%                                  6 and 7
%                  segparam(i,5) = turn-rate in degrees-per-second; this
%                                  number is positive for right turns and
%                                  negative for left turns; this
%                                  parameter is ignored for all segment
%                                  types other than 8 and 9
%                  segparam(i,6) = final heading in degrees; this
%                                  parameter is ignored for segment
%                                  types 6 and 8
%                  segparam(i,7) = turn-delta in degrees; this number is
%                                  for turn anticipation; specifically, the
%                                  turn generally must be stopped prior to
%                                  reaching the desired final heading since
%                                  there will be some overshoot as the
%                                  aircraft pitches and rolls back to level
%                                  flight; the turn-delta is the amount of
%                                  turn anticipation; it is found by trial-
%                                  and-error but can be ignored (i.e., set
%                                  to zero) if the exact value of the final
%                                  heading (i.e., the exact amount of the 
%                                  turn) is not considered critical; this
%                                  parameter is ignored for all segment
%                                  types other than 9
%                  segparam(i,8) = time step for i-th segment (seconds)
%
%   OUTPUTS
%       profile = flight profile
%          profile(i,1:3) = ENU path generated; 1=x, 2=y, 3=z
%          profile(i,4:6) = ENU velocity; 4 = x-velocity,
%                           5 = y-velocity, 6 = z-velocity 
%          profile(i,7:9) = NOT USED
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

%	M. & S. Braasch 02-2005
%	Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%
global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha
global thtl_final el_final ail_final rdr_final
x = [0 0 0 0 0 0 0 0 0 0 0 0 0]';
xd = [0 0 0 0 0 0 0 0 0 0 0 0 0]';
time = 0;
xcg = 0.3;

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
   deltat = segparam(i,8);
   
   if i == 1,
      pos = initpos;
      vel = initvel;
      dcm = initdcm;
   else,
      k = size(profile,1);
      pos = profile(k,1:3);
      vel = profile(k,4:6);
      dcm(1,1:3) = profile(k,10:12);
      dcm(2,1:3) = profile(k,13:15);
      dcm(3,1:3) = profile(k,16:18);
   end,

   if segparam(i,1) == 5,
      duration = segparam(i,2);
      v_final = segparam(i,3);
      heading = segparam(i,6)*pi/180;
      profilez = pro_f16strai(pos,vel,heading,duration,v_final,deltat);
   elseif segparam(i,1) == 6,
      gamma = segparam(i,4)*pi/180;  % global variable
      profilez = pro_f16pitch(pos,vel,dcm,deltat);
   elseif segparam(i,1) == 7,
      gamma = segparam(i,4)*pi/180;  % global variable
      duration = segparam(i,2);
      heading = segparam(i,6)*pi/180;
      profilez = pro_f16climb(pos,vel,heading,duration,deltat);
   elseif segparam(i,1) == 8,
      turnrate = segparam(i,5)*pi/180;  % global variable
      profilez = pro_f16pitroll(pos,vel,dcm,deltat);
   elseif segparam(i,1) == 9,
      turnrate = segparam(i,5)*pi/180;  % global variable
      hdg_final = segparam(i,6)*pi/180;
        if hdg_final < 0, hdg_final = hdg_final + 2*pi; end
      turndelta = segparam(i,7)*pi/180;
      profilez = pro_f16turn(pos,vel,dcm,hdg_final,turndelta,deltat);
   end
   
   if i == 1,
      mm = size(profilez(:,18),1);
      profilez(:,19) = [0:deltat:deltat*(mm-1)]';
      profile = profilez;
   else,
      mm = size(profilez(:,18),1);
      nn = size(profile,1);
      oldtime = profile(nn,19);
      profilez(1:mm,19) = oldtime + deltat*(1:mm)';
      profile = [profile; profilez];
   end
   
end

