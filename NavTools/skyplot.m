function skyplot(svxyzmat,svid,usrxyz,figflg,idflg)
%SKYPLOT        Generate satellite azimuth/elevation angle plot
%
%       SKYPLOT(svxyzmat,svid,usrxyz) or SKYPLOT(svenumat,svid) or
%          SKYPLOT(svxyzmat,svid,usrxyz,figflg,idflg)
%
%   INPUTS
%       svxyzmat = matrix of satellite positions in cartesian ECEF coordinates
%                  svxyzmat(i,1:3) = x,y,z coordinates for satellite i
%       svid = vector of satellite identification numbers corresponding to the
%              satellite locations given in SVXYZMAT
%       usrxyz = user position in cartesian ECEF coordinates
%
%       svenumat = matrix of satellite positions in East-North-Up coordinates
%                  relative to the user-determined skyplot origin (usually user location)
%       
%       figflg = 1 if the Matlab function 'close' is to be executed before plotting;
%                0 if the 'close' function is NOT to be executed before plotting;
%                Default value is 1
%       idflg = 1 if satellite identification numbers are desired on the plot,
%               0 otherwise.  Default is 1

%	References: 
%                   Understanding GPS: Principles and Applications,
%	            Elliott D. Kaplan, Editor, Artech House Publishers,
%	            Boston, 1996.
%
%	M. & S. Braasch 12-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.
%

if nargin<5, idflg=1; end
if nargin<4, figflg=1; end
flag=1;
if nargin<3, flag=0; end
if nargin<2, error('insufficient number of input arguments'),end

rad2deg=180/pi;
[N,M] = size(svxyzmat);
if flag,
   for i=1:N,
      matsvenu(i,:)=xyz2enu(svxyzmat(i,:),usrxyz)';
   end
else,
   matsvenu=svxyzmat;
end
for i=1:N,
  a=(pi/2)-atan2(matsvenu(i,1),matsvenu(i,2));
  r=0.01*abs( rad2deg*atan2(matsvenu(i,3),norm(matsvenu(i,1:2))) - 90 );
  svx(i)=r*cos(a); svy(i)=r*sin(a);
end
hlinex=[-0.9 0.9]; hliney=[0 0]; vlinex=[0 0]; vliney=[-0.9 0.9];
arg=[0:100]*2*pi/100;
if figflg,
   close
end
plot(0.9*sin(arg),0.9*cos(arg),'-',...
     0.8*sin(arg),0.8*cos(arg),'-',...
     0.7*sin(arg),0.7*cos(arg),'-',...
     0.6*sin(arg),0.6*cos(arg),'-',...
     0.5*sin(arg),0.5*cos(arg),'-',...
     0.4*sin(arg),0.4*cos(arg),'-',...
     0.3*sin(arg),0.3*cos(arg),'-',...
     0.2*sin(arg),0.2*cos(arg),'-',...
     0.1*sin(arg),0.1*cos(arg),'-',...
     hlinex,hliney,'-',vlinex,vliney,'-',...
     svx,svy,'*r')
axis('square')
axis('off')
title('SATELLITE SKYPLOT')
text(1.0,0,'EAST')
text(-1.3,0,'WEST')
text(-0.1,1,'NORTH')
text(-0.1,-1,'SOUTH')
if idflg,
   for i=1:N,
     text(svx(i)+0.05,svy(i),num2str(svid(i)))
   end
end
