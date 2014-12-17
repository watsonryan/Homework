%  sndemo11.m       Example of GEO satellite position generation and plotting
%
%       The location (50 degrees North, 5 degrees East) has been chosen
%       in the middle of Western Europe to give an indication of the
%       geo coverage provided by the WAAS and EGNOS
%
%       These are the geos:
%                120:  INMARSAT, Atlantic Ocean Region - East (AOR-E)
%                122:  INMARSAT, Atlantic Ocean Region - West (AOR-W)
%                124:  ARTEMIS
%                131:  INMARSAT, Indian Ocean Region (IOR)
clear all
close all
t = 0;
usrllh = [50*pi/180 5*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgeo
[svxyzmat,svid] = gensv(usrxyz,t,5);
skyplot(svxyzmat,svid,usrxyz)
text(0.6,0.8,'GEOs only')
text(-1.5,-1,'User at Lat 50 deg, Lon 5 deg')