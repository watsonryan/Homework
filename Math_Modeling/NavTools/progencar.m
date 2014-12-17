function profile = progencar(initpos,initvel,initdcm,segparam)
%PROGENCAR  	Land vehicle path (trajectory) generator.  
%               Path consists of straight segments joined by 
%               constant-radius turns.
%       
%	profile = progencar(initpos,initvel,initdcm,segparam)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (East and North 
%                 cartesian coordinates) (m/s)
%                 Note that the magnitude of the linear velocity 
%                 is constant throughout the trajectory
%       initdcm = initial direction cosine matrix (nav-to-body)
%                 for vehicle attitude (3x3 matrix)
%       segparam = segment and turn parameters;  N x 6 matrix
%                  where ...
%                  segparam(i,1) = segment-type identifier
%                       11 = level, straight constant-acceleration
%                           (including constant-velocity as a special case)
%                       12 = level constant-radius turn
%                  segparam(i,2) = duration (in seconds) of i-th 
%                                  straight segment (set to zero if
%                                  no straight portion is desired in
%                                  this segment).  This parameter is
%                                  ignored for segment type 12
%                  segparam(i,3) = final velocity in meters-per-second; this
%                                  parameter is ignored for segment type 12
%                  segparam(i,4) = direction and amount of turn (degrees).
%                                  Positive is a right turn.  This
%                                  parameter is ignored for segment type 11
%                  segparam(i,5) = turn-rate in degrees-per-second; this
%                                  parameter is ignored for segment type 11
%                  segparam(i,6) = time-step in seconds for i-th segment
%
%   OUTPUTS
%       profile = vehicle path profile
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
%   NOTE:  Only horizontal paths are supported.  The height specified 
%          by INITPOS will be held constant throughout the generated path.
%

%	M. & S. Braasch 03-2005
%	Copyright (c) 2005 by GPSoft
%	All Rights Reserved.
%

if nargin<4,error('insufficient number of input arguments'),end

vel(3)=0;

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
   deltat = segparam(i,6);
   
   if i == 1,
      pos = initpos;
      vel = initvel;
      dcmnb = initdcm;
   else
      k = size(profile,1);
      pos = profile(k,1:3);
      vel = profile(k,4:6);
      dcmnb(1,1:3) = profile(k,10:12);
      dcmnb(2,1:3) = profile(k,13:15);
      dcmnb(3,1:3) = profile(k,16:18);
   end,

   if segparam(i,1) == 11,
      duration = segparam(i,2);
      v_final = segparam(i,3);
      profilez = pro_car_strai(pos,vel,dcmnb,duration,v_final,deltat);
   elseif segparam(i,1) == 12,
      turnamt = segparam(i,4)*pi/180;
      turnrate = segparam(i,5)*pi/180;
      profilez = pro_car_turn(pos,vel,dcmnb,turnamt,turnrate,deltat);
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


