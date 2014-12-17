%  sndemo06.m
%
%     Baseline parity-space based 
%     fault detection example:  No failures
%
clear all
close all
%    
mpmat=mpgen(24,3600);
usrllh = [0*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
i=0;  initpos = [0 0 0 0];
randn('state',309874);
bar1 = waitbar(0,'Finding Residuals...  ');
for t = 1000:5:2800,
    i=i+1;
    [svxyzmat,svid] = gensv(usrxyz,t);
    [prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmat);
    estusr = olspos(prvec,svxyzmat);
    parvec = parityvec(prvec,svxyzmat,estusr);
    p(i) = norm(parvec);
    numsv(i) = max(size(prvec));
    svident(i) = svid(1);
    time(i) = t;
    waitbar(i/360)
end
close(bar1);

plot(time,p)
title('Baseline Fault Detection Scenario:  Normal Operation')
ylabel('magnitude of parity vector')
xlabel('GPS time of week in seconds')
