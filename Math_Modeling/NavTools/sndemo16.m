%  sndemo16.m   Example of GPS+GEO+Galileo satellite position generation and plotting
%
%   Extension of Demo 11
%
clear all
close all
t = 0;
usrllh = [50*pi/180 5*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
loadgeo
loadgalileo
[svxyzmat,svid] = gensv(usrxyz,t,5);
skyplot(svxyzmat,svid,usrxyz)
text(0.6,0.8,'GPS+Galileo+GEOs')
text(-1.5,-1,'User at Lat 50 deg, Lon 5 deg')