%  sndemo18.m  GPS+GEO+Galileo (Assumes common time base for all)
%              This is an extension of SNDEMO17.  We now assume
%              that WAAS and EGNOS allow us to reduce the effect
%              of the ionosphere by 75% (note the change of parameters
%              in the call to GENRNG)
clear all
close all
%
disp('    ')
disp('         Generating data  -   Please be patient...')
disp('   ')
mpmat=mpgen(230,3600,1,54321);
usrllh = [50*pi/180 5*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
loadgeo
loadgalileo
i=0;
bar1 = waitbar(0,'Calculating GNSS Position...  ');
randn('state',74347098);
for t = 41000:10:42800,
    i=i+1;
    [svxyzmat,svid] = gensv(usrxyz,t,5);
    [prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 0.25],[],mpmat);
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
title('WAAS/EGNOS Position Error  -  IONO Reduction Only')
ylabel('north error (m)')
xlabel('east error (m)')
