function profile = pro_car_strai(initpos,initvel,initdcm,duration,v_final,deltat)
%PRO_CAR_STRAI       Automobile profile sub-generator for a straight-and-level,
%                    constant acceleration or constant velocity path segment.
%       
%   profile = pro_car_strai(initpos,initvel,initdcm,duration,v_final,deltat);
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (3 ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (3 ENU cartesian 
%                 coordinates) (m/s)
%       initdcm = initial direction cosine matrix for
%                 vehicle attitude (navigation to body
%                 frame) (3x3 matrix)
%       duration = length of flight segment (seconds)
%       v_final = total final velocity (i.e., speed) in meters-per-second
%       deltat = time increment in seconds
%
%   OUTPUTS
%       profile = vehicle path profile
%          profile(i,1:3) = ENU path generated; 1=x, 2=y, 3=z;
%                           units are meters
%          profile(i,4:6) = ENU velocity; 4 = x-velocity,
%                           5 = y-velocity, 6 = z-velocity;
%                           units are meters/second
%          profile(i,7:9) = NOT USED (format retained for compatibility
%                           with other toolbox functions
%          profile(i,10:18) = elements of the direction cosine matrix
%                            (DCM) for vehicle attitude; 10 = DCM(1,1),
%                            11 = DCM(1,2), 12 = DCM(1,3),
%                            13 = DCM(2,1), et cetera
%

%	M. & S. Braasch 03-2005
%	Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<6,error('insufficient number of input arguments'),end

dcmbn=initdcm';
eulvect=dcm2eulr(dcmbn);
%%phi=eulvect(1);
%%theta=eulvect(2);
psi=eulvect(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = deltat:deltat:duration; t=t';
    tmp = ones(max(size(t)),1);
    
    gnd_trk_ang = psi;
    v_init = norm(initvel);
    accel = (v_final - v_init)/duration;

    profile(:,7) = NaN*tmp;
    profile(:,8) = NaN*tmp;
    profile(:,9) = NaN*tmp;

    profile(:,4) = v_init*sin(gnd_trk_ang)*tmp + accel*sin(gnd_trk_ang)*t;
    profile(:,5) = v_init*cos(gnd_trk_ang)*tmp + accel*cos(gnd_trk_ang)*t;
    profile(:,6) = 0;

    profile(:,1) = initpos(1)*tmp + v_init*sin(gnd_trk_ang)*t + 0.5*accel*sin(gnd_trk_ang)*t.*t;
    profile(:,2) = initpos(2)*tmp + v_init*cos(gnd_trk_ang)*t + 0.5*accel*cos(gnd_trk_ang)*t.*t;
    profile(:,3) = initpos(3)*tmp;
    
    profile(:,10) = initdcm(1,1)*tmp;
    profile(:,11) = initdcm(1,2)*tmp;
    profile(:,12) = initdcm(1,3)*tmp;
    profile(:,13) = initdcm(2,1)*tmp;
    profile(:,14) = initdcm(2,2)*tmp;
    profile(:,15) = initdcm(2,3)*tmp;
    profile(:,16) = initdcm(3,1)*tmp;
    profile(:,17) = initdcm(3,2)*tmp;
    profile(:,18) = initdcm(3,3)*tmp;

    