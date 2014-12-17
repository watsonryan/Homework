%  sndemo02.m         Short example of satellite motion through sky
clear all
close all
%    
usrllh = [0*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps

[svxyzmat,svid] = gensv(usrxyz,0,0);		%time = 0 case places svid#
pause(0.1)
skyplot(svxyzmat,svid,usrxyz,0,1)
hold on


for t = 200:200:3400,
  [svxyzmat,svid] = gensv(usrxyz,t,0);
  pause(0.1)
  skyplot(svxyzmat,svid,usrxyz,0,0)
  hold on
end

[svxyzmat,svid] = gensv(usrxyz,3600,0);		%time = 3600 case places svid#
pause(0.1)
skyplot(svxyzmat,svid,usrxyz,0,1)
hold on

text(0.6,0.9,'Satellite motion depicted')
text(0.6,0.8,'over 1 hour')
text(-1.5,-0.9,'User at Lat 0 deg, Lon 0 deg')