%  sndemo15.m     GPS + Galileo 
%             (incorporates possible GPS/Galileo system time offset)
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
loadgalileo
toffset = 300;  % time offset in meters
i=0;
bar1 = waitbar(0,'Generating GNSS Position...  ');
randn('state',74347098);
for t = 41000:10:42800,
    i=i+1;
    [svxyzmat,svid] = gensv(usrxyz,t,5);
    [prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmat);
    for j = 1:max(size(svid)),
        if svid(j) > 49, prvec(j) = prvec(j) + toffset; end
    end
    estusr = olsposgg(prvec,svxyzmat,svid);
    enuerr(i,:) = ( xyz2enu(estusr(1:3),usrxyz) )';
    rxterr(i) = estusr(4);  % true clk bias is zero
    gpsgalterr(i) = estusr(5) - toffset;
    waitbar(i/180)
end
close(bar1);
close
plot(enuerr(:,1),enuerr(:,2),'*')
axis('square')
axis('equal')
axis([-10 10 -10 10])
grid
title('GPS + Galileo Position Error  -  Half-Hour Scenario')
ylabel('north error (m)')
xlabel('east error (m)')
text(-4,7,'Non-Zero GPS/Galileo Time Offset Assumed')
