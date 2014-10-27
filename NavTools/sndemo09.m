%  sndemo09.m   Compute DOP's and compare with
%               Monte Carlo trials
clear all
close all
%    
time = 0;
usrllh = [0*pi/180 0*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
[svxyzmat,svidvec] = gensv(usrxyz,time,15);  % The satellite positions will be held
                                             % fixed for these Monte Carlo
                                             % trials
bar1 = waitbar(0,'Running Monte Carlo Trials...  ');
randn('state',7912345);
for i = 1:1000,
    [prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svidvec,time,[1 0 0 0 0]);
    estusr = olspos(prvec,svxyzmat);
    err(i,1:3) = (xyz2enu(estusr(1:3),usrxyz))';  %estusr(1:3) - usrxyz(1:3);
    terr(i) = estusr(4);  % true clk bias is zero
    waitbar(i/1000)
 end
close(bar1);
for j = 1:length(svidvec),
    svenumat(j,:) = (xyz2enu(svxyzmat(j,:),usrxyz))';
end
dopvec = dops(svenumat,[0 0 0]);
subplot(311)
plot(1:1000,err(:,1))
axis([1 1000 -4 4])
ylabel('east error [m]')
title('Monte Carlo Results: Fixed Satellite Locations')
subplot(312)
plot(1:1000,err(:,2))
axis([1 1000 -4 4])
ylabel('north error [m]')
subplot(313)
plot(1:1000,err(:,3))
axis([1 1000 -4 4])
ylabel('up error [m]')
xlabel('Monte Carlo trial number')

dopvec(1:3)
stdxerr = std(err(:,1))
stdyerr = std(err(:,2))
stdzerr = std(err(:,3))
disp('ANS above is listing the first three values of DOPVEC')
disp('which are the theoretical XDOP, YDOP and ZDOP values.')
disp('Compare these to the estimated values from the Monte')
disp('Carlo trials: stdxerr, stdyerr and stdzerr.')
