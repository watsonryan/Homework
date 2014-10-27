function profile = prostrai(initpos,initvel,initacc,initdcm,duration,deltat)
%PROSTRAI        Flight profile sub-generator for a straight,
%               constant velocity or constant acceleration
%               flight segment (no turn).  Local-level (i.e., East-
%               North-Up coordinates) version suitable for 
%               short-distance, short-duration flights.
%       
%	profile = prostrai(initpos,initvel,initacc,initdcm,duration,deltat)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (3 ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (3 ENU cartesian 
%                 coordinates) (m/s)
%       initacc = initial acceleration vector (m/s^2)
%       initdcm = initial direction cosine matrix for
%                 vehicle attitude (3x3 matrix)
%       duration = length of flight segment (seconds)
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
%

%	M. & S. Braasch 12-97
%	Copyright (c) 1997 by GPSoft
%	All Rights Reserved.
%

if nargin<5,deltat=1;end
if nargin<4,error('insufficient number of input arguments'),end

    t = deltat:deltat:duration; t=t';
    tmp = ones(max(size(t)),1);

    profile(:,7) = initacc(1)*tmp;
    profile(:,8) = initacc(2)*tmp;
    profile(:,9) = initacc(3)*tmp;

    profile(:,4) = initvel(1)*tmp + initacc(1)*t;
    profile(:,5) = initvel(2)*tmp + initacc(2)*t;
    profile(:,6) = initvel(3)*tmp + initacc(3)*t;

    profile(:,1) = initpos(1)*tmp + initvel(1)*t + 0.5*initacc(1)*t.*t;
    profile(:,2) = initpos(2)*tmp + initvel(2)*t + 0.5*initacc(2)*t.*t;
    profile(:,3) = initpos(3)*tmp + initvel(3)*t + 0.5*initacc(3)*t.*t;

    profile(:,10) = initdcm(1,1)*tmp;
    profile(:,11) = initdcm(1,2)*tmp;
    profile(:,12) = initdcm(1,3)*tmp;
    profile(:,13) = initdcm(2,1)*tmp;
    profile(:,14) = initdcm(2,2)*tmp;
    profile(:,15) = initdcm(2,3)*tmp;
    profile(:,16) = initdcm(3,1)*tmp;
    profile(:,17) = initdcm(3,2)*tmp;
    profile(:,18) = initdcm(3,3)*tmp;
