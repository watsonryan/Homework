%  sndemo38.m   Demo of Satellite Visibility using Rinex Ephemeris data
%
clear all
close all
%
mthnxyz = [901906.698 -5726170.573 2651517.227];   % Marathon, FL  CORS station
filename = 'day101s.03n';
daynum = 0;
tstart = 18*3600;  tstop = 20*3600;  tinc = 60;
plotflg = 2;  maskang = 5;
vismat = satvis2(mthnxyz,tstart,tstop,tinc,daynum,filename,plotflg,maskang);
