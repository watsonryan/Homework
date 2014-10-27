function profile = pro_car_turn(initpos,initvel,initdcm,turnamt,turnrate,deltat)
%PRO_CAR_TURN       Automobile profile sub-generator for a constant-radius
%                   turn segment.
%       
%   profile = pro_car_turn(initpos,initvel,initdcm,turnamt,turnrate,deltat)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (3 ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (3 ENU cartesian 
%                 coordinates) (m/s)
%       initdcm = initial direction cosine matrix for
%                 vehicle attitude (navigation to body
%                 frame) (3x3 matrix)
%       turnamt = direction and amount of turn (radians).
%                 Positive is a right turn.
%       turnrate = turn-rate in radians-per-second
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

dcm = initdcm;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % constant-radius turn
    vel = initvel;
    r = norm(initvel)/turnrate;
    numstps=round((abs(turnamt)/turnrate)/deltat);
        if numstps<2, 
           fprintf(1,'turn cannot be simulated !\n',j)
           error('specify smaller time step or smaller turn rate')
        end
    anginc=(turnamt)/numstps;
    velang=atan2(initvel(2),initvel(1));
    beta=velang-(pi/2)*sign(turnamt);
    turnorg(1)=initpos(1)+r*cos(beta);
    turnorg(2)=initpos(2)+r*sin(beta);
    cumang=0;
 
    rotmat=[cos(anginc) sin(anginc) 0;
            -sin(anginc) cos(anginc) 0;
               0             0       1];

    h = waitbar(0,'Generating turn segment');
    for i = 1:numstps,
        cumang = cumang + anginc;
        path(1)=turnorg(1)+r*cos(beta+pi-cumang);
        path(2)=turnorg(2)+r*sin(beta+pi-cumang);
        path(3)=initpos(3);

        vel = (rotmat*vel')';
        dcm = rotmat*dcm;
        
        profile(i,7) = NaN;
        profile(i,8) = NaN;
        profile(i,9) = NaN;

        profile(i,4) = vel(1);
        profile(i,5) = vel(2);
        profile(i,6) = 0;

        profile(i,1) = path(1);
        profile(i,2) = path(2);
        profile(i,3) = path(3);
    
        profile(i,10) = dcm(1,1);
        profile(i,11) = dcm(1,2);
        profile(i,12) = dcm(1,3);
        profile(i,13) = dcm(2,1);
        profile(i,14) = dcm(2,2);
        profile(i,15) = dcm(2,3);
        profile(i,16) = dcm(3,1);
        profile(i,17) = dcm(3,2);
        profile(i,18) = dcm(3,3);
        
        waitbar(i/numstps,h)

    end,
    close(h)
