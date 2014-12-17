function path = pathgen(initpos,initvel,segparam,deltat)
%PATHGEN	Vehicle path (trajectory) generator.  
%               Path consists of straight segments joined by 
%               constant-radius turns.
%       
%	path = pathgen(initpos,initvel,segparam,deltat)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (East and North 
%                 cartesian coordinates) (m/s)
%                 Note that the magnitude of the linear velocity 
%                 is constant throughout the trajectory
%       segparam = segment and turn paramters
%                  segparam(i,1) = duration (in seconds) of i-th 
%                                  straight segement
%                  segparam(i,2) = direction and amount of i-th turn 
%                                  (set to zero if no turn is desired) 
%                                  (degrees) Positive is a left turn.
%                  segparam(i,3) = centripetal acceleration of i-th 
%                                  turn (in m/s^2)
%       deltat = time increment in seconds
%
%   OUTPUTS
%	path = ENU path generated
%
%   NOTE:  Only horizontal paths are supported.  The height specified 
%          by INITPOS will be held constant throughout the generated path.

%	M. & S. Braasch 12-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.
%

if nargin<4,deltat=1;end
if nargin<3,error('insufficient number of input arguments'),end

[m,n]=size(initvel); if m>n, vel=initvel'; else, vel=initvel; end
vel(3)=0;
[m,n]=size(initpos); if m>n, initpos=initpos'; end
[m,n]=size(segparam);
i=1;

path(i,:)=initpos;
for j = 1:m,
    % straight segment
    for t = deltat:deltat:segparam(j,1),
        i = i + 1;
        path(i,:) = path(i-1,:) + vel*deltat; 
    end

    % constant-radius turn
    if segparam(j,2) ~= 0,
        r = (norm(vel))^2/segparam(j,3);
        omega = norm(vel)/r;
        numstps=round(((abs(segparam(j,2))*pi/180)/omega)/deltat);
        if numstps<2, 
           fprintf(1,'turn number %i cannot be simulated !\n',j)
           error('specify smaller time step or smaller centripetal acceleration')
        end
        anginc=(segparam(j,2)*pi/180)/numstps;
        velang=atan2(vel(2),vel(1));
        beta=velang+(pi/2)*sign(segparam(j,2));
        turnorg(1)=path(i,1)+r*cos(beta);
        turnorg(2)=path(i,2)+r*sin(beta);
        cumang=0;
        for kk = 1:numstps,
            i=i+1;
            cumang = cumang + anginc;
            path(i,1)=turnorg(1)+r*cos(beta-pi+cumang);
            path(i,2)=turnorg(2)+r*sin(beta-pi+cumang);
            path(i,3)=path(i-1,3);
        end,
        rotmat=[cos(-cumang) sin(-cumang) 0;
                -sin(-cumang) cos(-cumang) 0;
                   0             0       1];
        vel=(rotmat*vel')';
    end
end
