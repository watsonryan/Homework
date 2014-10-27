%  sndemo20.m    Satellite Visibility
%
clear all
close all
%
disp('    ')
disp('         Generating data  -   Please be patient...')
disp('   ')
usrllh = [50*pi/180 5*pi/180 0];
usrxyz = llh2xyz(usrllh);
tstart = 0;  tstop = 86000;  tinc = 100;  sysflg = [1 3 4];
plotflg = 1;  maskang = 5;
vismat = satvis(usrxyz,tstart,tstop,tinc,sysflg,plotflg,maskang);
text(10,10,'GPS+Galileo+GEOs')
