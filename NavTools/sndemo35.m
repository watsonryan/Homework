%  sndemo35.m       Static User   OLS positioning
%              SAME AS SNDEMO05.M BUT THIS TIME USE YUMA-FORMATTED
%              ALMANAC DATA
clear all
close all
%    
mpmat=mpgen(32,3600,1,54321);    % Must generate for 32 satellites
%                                % since satellite id numbers range
%                                % from 1 to 32 in the actual almanac
usrllh = [0*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadyuma('yuma1.txt')           % Load the Yuma almanac file
i=0;
bar1 = waitbar(0,'Calculating Position...  ');
randn('state',74347098);
for t = 41000:10:42800,
    i=i+1;
    [svxyzmat,svid] = gensvalm(usrxyz,t,5);   % Use the almanac-specific
    %                                         % function for SV positioning
    [prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmat);
    estusr = olspos(prvec,svxyzmat);
    enuerr(i,:) = ( xyz2enu(estusr(1:3),usrxyz) )';
    terr(i) = estusr(4);  % true clk bias is zero
    waitbar(i/180)
end
close(bar1); 
close
plot(enuerr(:,1),enuerr(:,2),'*')
axis('square')
axis('equal')
axis([-10 10 -10 10])
grid
title('GPS Positioning Error  -  Static User  -  Yuma-Format Almanac Data')
ylabel('north error (m)')
xlabel('east error (m)')
