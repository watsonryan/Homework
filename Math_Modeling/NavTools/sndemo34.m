%  sndemo34.m       Example of satellite position generation and plotting
%              SAME AS SNDEMO01.M BUT THIS TIME USE YUMA-FORMATTED
%              ALMANAC DATA
clear all
close all
%    
t = 0;
usrllh = [0*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadyuma('yuma1.txt')                         % Load the Yuma almanac file
[svxyzmat,svid] = gensvalm(usrxyz,t,0);    % Use the almanac-specific
%                                          % function for sv positioning
skyplot(svxyzmat,svid,usrxyz)
text(-1.6,-1.1,'Yuma-Formatted Almanac Data Utilized')