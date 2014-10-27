%  sndemo26.m    Skyplot
%
clear all
close all
%
usrllh = [0.1*pi/180 0.1*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
t = 1000;
[svmatusr,svidusr] = gensv(usrxyz,t,5);
skyplot(svmatusr,svidusr,usrxyz)
