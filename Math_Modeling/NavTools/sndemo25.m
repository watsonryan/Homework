%  sndemo25.m    Satellite Visibility
%
clear all
close all
%
usrllh = [0.1*pi/180 0.1*pi/180 0];
usrxyz = llh2xyz(usrllh);
tstart = 1000;  tstop = 2800;  tinc = 10;  sysflg = 1;
plotflg = 2;  maskang = 5;
vismat = satvis(usrxyz,tstart,tstop,tinc,sysflg,plotflg,maskang);
