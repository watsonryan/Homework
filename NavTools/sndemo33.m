%  sndemo33.m    Multipath Analysis
%
clear all
close all
%
mpmat=mpgen(24,3600,1,33333);
refllh = [0*pi/180 0*pi/180 0];
refxyz = llh2xyz(refllh);

loadgps
i=0;
bar1 = waitbar(0,'Generating Ranges...  ');
randn('state',74347098);
for t = 41001:1:41600,
    i=i+1;
    [svmatref,svidref] = gensv(refxyz,t,0);
    [prref,adrref] = genrng(1,refxyz,svmatref,svidref,t,...
                            [1 1 0 1 1],[],mpmat);
    pr8(i) = prref(find(svidref==8));
    adr8(i) = adrref(find(svidref==8));
	 waitbar(i/600)    
end
close(bar1); 
resi = pr8 - adr8;
resi = resi - mean(resi);
plot(resi)
title('Multipath, Thermal Noise and Iono Divergence')
ylabel('residuals in meters')
xlabel('run time in seconds')
