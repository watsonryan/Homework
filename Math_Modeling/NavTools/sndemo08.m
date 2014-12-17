%  sndemo08.m
%
%     PARITYVEC example:  1 m/s ramp error injected onto SV 3 at t = 1750 
%                         (750 seconds after start)
%
clear all
close all
%    
mpmat=mpgen(24,3600);
usrllh = [0*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
i=0;  initpos = [0 0 0 0];  ramperr = 0;
randn('state',309874);
bar1 = waitbar(0,'Finding Residuals...  ');
for t = 1000:5:2800,
    i=i+1;
    [svxyzmat,svid] = gensv(usrxyz,t);
    [prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmat);
    if t > 1750,
       ramperr = 1*(t-1750);
       prvec(1) = prvec(1) + ramperr;
    end
    estusr = olspos(prvec,svxyzmat);
    sv(i) = svid(1);
    parvec = parityvec(prvec,svxyzmat,estusr);
    p(i) = norm(parvec);
    time(i) = t;
    terr(i) = estusr(4);  % true clk bias is zero
    waitbar(i/360)
 end
close(bar1);
plot(time,p)
title('Fault Detection Example: Ramp error on SV 3 Injected at t = 1750')
ylabel('magnitude of parity vector')
xlabel('GPS time of week in seconds')
