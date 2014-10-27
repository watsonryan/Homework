%  sndemo01.m       Example of satellite position generation and plotting
clear all
close all
%    
t = 39600;
usrllh = [40*pi/180 -90*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
[svxyzmat,svid] = gensv(usrxyz,t,5);
skyplot(svxyzmat,svid,usrxyz)
text(0.5,0.9,'Simulated GPS Constellation')
text(0.5,0.8,'User at Lat 40 deg, Lon -90 deg')
