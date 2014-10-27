function [profile,errflg] = ...
   proturn(initpos,initvel,initacc,inidcmnb,turnamt,deltat)
%PROTURN        Flight profile sub-generator for a constant altitude,
%               constant acceleration turn.  Local-level (i.e., East-
%               North-Up coordinates) version suitable for 
%               short-distance, short-duration flights.
%       
%  [profile,errflg] = proturn(initpos,initvel,initacc,inidcmnb,turnamt,deltat)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (3 ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector; z component must
%                 zero since this is a level turn (3 ENU cartesian 
%                 coordinates) (m/s)
%       initacc = initial acceleration vector (this will be zeros
%                 since this is a constant velocity turn; it is
%                 included here for consistency with other
%                 progenX functions (m/s^2)
%       inidcmnb = initial direction cosine matrix for
%                 vehicle attitude (navigation to body
%                 frame) (3x3 matrix)
%       turnamt = amount of turn (degrees)
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
%               = 1 if turn cannot be simulated and the user should
%                 specify a smaller time step or smaller centripetal
%                 acceleration
%               = 2 if bank angle is less than 1 degree
%
%   NOTES
%       A constant bank angle (as given by INIDCMNB) is used 
%       throughout the turn.  Must use PROROLL to provide a
%       transition between turn and straight-and-level 
%       (both before and after turn).

%	M. & S. Braasch 12-97
%	Copyright (c) 1997 by GPSoft
%	All Rights Reserved.
%

if nargin<6,deltat=1;end
if nargin<5,error('insufficient number of input arguments'),end
[m,n]=size(initpos); if m>n, pos=initpos'; else, pos=initpos; end
[m,n]=size(initvel); if m>n, vel=initvel'; else, vel=initvel; end
vel(3) = 0;
errflg = 0;

dcmbn=inidcmnb';
eulvect=dcm2eulr(dcmbn);
phi=eulvect(1); theta=eulvect(2); psi=eulvect(3);
if abs(phi) < 1*pi/180,
   error('bank angle less than 1 degree:  too small for a turn!')
end

cntrpacc = abs(9.81*tan(phi));

profile(1,1:3) = pos;

        r = (norm(vel))^2/cntrpacc;
        omega = norm(vel)/r;
        numstps=round(((turnamt*pi/180)/omega)/deltat);
        if numstps<2, 
           errflg = 1;
           fprintf(1,'turn cannot be simulated !\n')
           error('specify smaller time step or smaller bank angle')
        end
        anginc=(-1)*sign(phi)*(turnamt*pi/180)/numstps;
        velang=atan2(vel(2),vel(1));
        beta=velang+(pi/2)*(-1)*sign(phi);
        turnorg(1)=pos(1,1)+r*cos(beta);
        turnorg(2)=pos(1,2)+r*sin(beta);
        cumang=0; i=0;
        h = waitbar(0,'Generating turn profile');
        for kk = 1:numstps,
            i=i+1;
            cumang = cumang + anginc;
            profile(i,1)=turnorg(1)+r*cos(beta-pi+cumang);
            profile(i,2)=turnorg(2)+r*sin(beta-pi+cumang);
            profile(i,3)=pos(3);
            
            rotmat=[cos(-anginc) sin(-anginc) 0;
                   -sin(-anginc) cos(-anginc) 0;
                        0             0       1];    
            vel=(rotmat*vel')';
            profile(i,4:6) = vel;

            accdvect= turnorg(1:2) - profile(i,1:2);
            naccdvec=accdvect/norm(accdvect);
            profile(i,7:8) = cntrpacc*naccdvec;
            profile(i,9) = 0;

            psi = psi - anginc;
            dcmnb = eulr2dcm([phi theta psi]);
            profile(i,10) = dcmnb(1,1);
            profile(i,11) = dcmnb(1,2);
            profile(i,12) = dcmnb(1,3);
            profile(i,13) = dcmnb(2,1);
            profile(i,14) = dcmnb(2,2);
            profile(i,15) = dcmnb(2,3);
            profile(i,16) = dcmnb(3,1);
            profile(i,17) = dcmnb(3,2);
            profile(i,18) = dcmnb(3,3);
            waitbar(kk/numstps,h)
        end,
        close(h)
        