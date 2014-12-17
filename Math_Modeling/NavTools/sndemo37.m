%  sndemo37.m   Demo of Satellite Visibility using Yuma almanac data
%
clear all
close all
%
usrxyz = [678454.119 -4893824.462 4020518.314];
filename = 'yuma1.txt';
daynum = 0;  currweek = 160;
tstart = 0;  tstop = 3*3600;  tinc = 60;
plotflg = 2;  maskang = 5;
vismat = satvis3(usrxyz,tstart,tstop,tinc,daynum,filename,plotflg,currweek,maskang);
