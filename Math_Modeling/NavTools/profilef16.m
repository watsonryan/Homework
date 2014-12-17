function profile = profilef16(initpos,initvtm,initpsi,flightseg)
%PROFILEF16		F-16 flight profile generator.  Feedback control
%               simulation is incorporated.  Local-level version
%               suitable for short-distance, short-duration flights.
%               Profile can consist of:
%
%		- altitude-hold segments with specified desired airspeed
%       - Climbs or descents to specified altitude at a given flight path
%       angle
%		- constant-altitude, coordinated turns to a specified heading (at a
%		given turn-rate)
%       - constant acceleration segment to simulate the path along the
%       runway prior to rotation (take-off)
%
%	profile = profilef16(initpos,initvtm,initpsi,flightseg)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (ENU cartesian coordinates) (meters)
%       initvtm = initial airspeed (m/s)
%       initpsi = initial yaw angle (radians)
%       flightseg = flight segment parameters;  N x 8 matrix
%                  where ...
%                  flightseg(i,1) = segment-type identifier
%                       5 = altitude hold (i.e., straight-and-level); a
%                       desired velocity is specified that the aircraft
%                       will acclerate/decelerate to if necessary
%                       7 = climb or descent to a desired altitude at a
%                       given flight path angle
%                       9 = level, coordinated turn at a given turn rate
%                      10 = acceleration along the runway prior to
%                          rotation; vehicle is level, and accelerates over
%                          the specified duration up to the specified final
%                          speed.  This segment must be used to simulate
%                          the acceleration of the vehicle along the runway
%                          since segment #5 assumes the vehicle is flying
%                          and thus cannot be used at slow speeds (i.e.,
%                          less than 50 m/s)
%                  flightseg(i,2) = duration (in seconds) of i-th segment;
%                                   this parameter is ignored for all
%                                   segment types other than 5
%                  flightseg(i,3) = final velocity in meters-per-second;
%                                   this parameter is ignored (i.e.,
%                                   previous value of velocity is held
%                                   constant) for turns and climbs
%                  flightseg(i,4) = flight path angle (i.e., climb angle) in
%                                  degrees; this parameter is ignored for
%                                  all segment type numbers other than 7
%                  flightseg(i,5) = turn-rate in degrees-per-second; this
%                                  number is positive for right turns and
%                                  negative for left turns; this
%                                  parameter is ignored for all segment
%                                  types other than 9
%                  flightseg(i,6) = final heading in degrees;
%                  flightseg(i,7) = final altitude in meters;
%                  flightseg(i,8) = time step for i-th segment (seconds)
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
%       For each segment specified, all eight parameters in FLIGHTSEG must
%       be specified.  Dummy arguments must have a value assigned.
%

%	M. & S. Braasch 12-2006; 09-2007
%	Copyright (c) 2006-2007 by GPSoft LLC
%	All Rights Reserved.
%
%   Revision History
%      September 2007:  Changed the initial conditional to check if the
%                       initial speed is less than 75 m/s.  If so then the
%                       trim routine is not run (the aircraft is assumed
%                       to be on the ground and thus pitch and roll are
%                       zero).  This conditional previously checked to see
%                       if the aircraft was stationary (speed = 0).
%
global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha
global thtl_final el_final ail_final rdr_final
global m2f f2m
x = [0 0 0 0 0 0 0 0 0 0 0 0 0]';
xd = [0 0 0 0 0 0 0 0 0 0 0 0 0]';
time = 0;
xcg = 0.3;
f2m = 0.3048;  m2f = 1/f2m;   % feet/meters conversions

if nargin<4,error('insufficient number of input arguments'),end

[m,n]=size(initpos); if m>n, initpos=initpos'; end
[m,n]=size(flightseg); nsegs = m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Determine initial trim conditions

if initvtm < 75,    % if the initial vehicle speed is less than rotation
%                   % speed then do not run the trim routine

    x = zeros(13,1);
    x(6) = initpsi;
    x(10) = initpos(2)*m2f;
    x(11) = initpos(1)*m2f;
    x(12) = initpos(3)*m2f;
    x_trim_start = x;
else
    rollrate=0; pitchrate=0; turnrate=0; coord=0; stab=0;
    thtl_init = 0;
    el_init = 0;
    ail_init = 0;
    rdr_init = 0;
    alpha_init = 0;
    beta_init = 0;

    % gamma is passed in via global variable
    gamma = 0;
    x = [initvtm*m2f 0 0 0 0 initpsi 0 0 0 0 0 0 0]';
    [init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
    thtl_start = thtl_final;
    el_start = el_final;
    ail_start = ail_final;
    rdr_start = rdr_final;
    x(6) = initpsi;
    x(10) = initpos(2)*m2f;
    x(11) = initpos(1)*m2f;
    x(12) = initpos(3)*m2f;
    x_trim_start = x;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nsegs,

   fprintf(1,' Creating profile for segment number %i \n',i)
   deltat = flightseg(i,8);
   
   if flightseg(i,1) == 10,
      v_final = flightseg(i,3);
      duration = flightseg(i,2);
      [profilez,thtl_end,el_end,ail_end,rdr_end] = ...
          profilef16start(v_final,deltat,duration);
        thtl_start = thtl_end;
        el_start = el_end;
        ail_start = ail_end;
        rdr_start = rdr_end;
   elseif flightseg(i,1) == 5,
      duration = flightseg(i,2);
      v_final = flightseg(i,3);
      heading = flightseg(i,6)*pi/180;
      altitude = flightseg(i,7);
      [profilez,thtl_end,el_end,ail_end,rdr_end] = ...
          profilef16strai(altitude,heading,duration,v_final,deltat,...
          thtl_start,el_start,ail_start,rdr_start);
        thtl_start = thtl_end;
        el_start = el_end;
        ail_start = ail_end;
        rdr_start = rdr_end;
   elseif flightseg(i,1) == 7,
      gamma = flightseg(i,4)*pi/180;  % climb angle is a global variable
      heading = flightseg(i,6)*pi/180;
      altitude = flightseg(i,7);
      [profilez,thtl_end,el_end,ail_end,rdr_end] = ...
          profilef16climb(altitude,heading,deltat,...
          thtl_start,el_start,ail_start,rdr_start);
        thtl_start = thtl_end;
        el_start = el_end;
        ail_start = ail_end;
        rdr_start = rdr_end;
   elseif flightseg(i,1) == 9,
      turnrate = flightseg(i,5)*pi/180;  % global variable
      altitude = flightseg(i,7);
      hdg_final = flightseg(i,6)*pi/180;
        if hdg_final < 0, hdg_final = hdg_final + 2*pi; end
      [profilez,thtl_end,el_end,ail_end,rdr_end] = ...
          profilef16turn(hdg_final,altitude,deltat,thtl_start,el_start,...
          ail_start,rdr_start);
        thtl_start = thtl_end;
        el_start = el_end;
        ail_start = ail_end;
        rdr_start = rdr_end;
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

