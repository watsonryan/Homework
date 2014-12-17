function skyplot3(filename,tow,usrxyz,currweek,maskangle,figflg,idflg)
%SKYPLOT3        Generate satellite azimuth/elevation angle plot
%               using Yuma almanac data
%
%       SKYPLOT3(filename,tow,usrxyz,currweek,maskangle,figflg,idflg)
%
%   INPUTS
%       filename = name of the almanac file; note the almanac must be in 
%                  Yuma format; the filename must be in quotes 
%                  (e.g., 'yuma19.txt')
%       tow = GPS time of week in seconds
%       usrxyz = user position in cartesian ECEF coordinates
%       currweek = GPS week number for the time specified (default is the week
%                  number specified for the first satellite listed in the almanac file)
%       maskangle = satellite elevation mask angle in degrees (default is 5 degrees)
%       figflg = 1 if the Matlab function 'close' is to be executed before plotting;
%                0 if the 'close' function is NOT to be executed before plotting;
%                Default value is 1
%       idflg = 1 if satellite identification numbers are desired on the plot,
%               0 otherwise.  Default is 1
%
%   See also:  SKYPLOT, SKYPLOT2

%
%	M. & S. Braasch 12-2002
%	Copyright (c) 2002 by GPSoft
%	All Rights Reserved.
%
weekflag = 0;
if nargin<7, idflg=1; end
if nargin<6, figflg=1; end
if nargin<5, weekflag = 1; end
if nargin<4, maskangle = 5; end
if nargin<3, error('insufficient number of input arguments'),end

global SVIDV MV OMGV RV INCLV TOEV HEALTHV ECCENV
global OMGDOTV ARGPERIV AF0V AF1V WEEKV
loadyuma(filename);

id = find(SVIDV);

if weekflag == 1,
    currweek = WEEKV(id(1));
end
rad2deg=180/pi;

j = 0;
for i=1:length(id),
      %N = id(i);
      N = i;
      svxyz = svposalm(RV(N),TOEV(N),MV(N),OMGV(N),INCLV(N),...
           tow,ECCENV(N),ARGPERIV(N),OMGDOTV(N),WEEKV(N),currweek);
        svenu = xyz2enu(svxyz,usrxyz);
        el = (180/pi)*atan(svenu(3)/norm(svenu(1:2)));
        if el >= maskangle,
            j = j + 1;
            matsvenu(j,:)=svenu';
            svid(j) = SVIDV(i);
        end
end

N = j;
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
