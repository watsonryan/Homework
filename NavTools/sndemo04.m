    %  sndemo04.m      Satellite Motion/Visibility:  24 hours near the north pole
clear all
close all
%    
usrllh = [89.9*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
for t = 0:2000:86000,
  [svxyzmat,svid] = gensv(usrxyz,t,0);
  pause(0.1)
  skyplot(svxyzmat,svid,usrxyz,0,0)
  hold on
end
hold off
text(0.6,0.9,'Satellite motion depicted')
text(0.6,0.8,'over 24 hours')
text(-1.5,-0.9,'User at Lat 89.9 deg, Lon 0 deg')