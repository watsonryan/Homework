%  sndemo36.m      DGPS (Range corrections)
%                  No Carrier-Smoothing
%           SAME AS SNDEMO14.M BUT THIS TIME
%           USE YUMA-FORMATTED ALMANAC DATA
clear all
close all
%
                                   % Since satellite id numbers
mpmatusr=mpgen(32,3600,1,22222);   % range from 1 to 32, must
mpmatref=mpgen(32,3600,1,33333);   % generate SA and multipath
%                                  % for 32 rather than 24 sv's
refllh = [0*pi/180 0*pi/180 0];
refxyz = llh2xyz(refllh);
usrllh = [0.1*pi/180 0.1*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadyuma('yuma1.txt')              % Specify which almanac file to use
i=0;
bar1 = waitbar(0,'Generating DGPS Position...  ');
for t = 1000:10:2800,
    i=i+1;
%%%Reference Station
    clear svxyzmat svid prvec adrvec
    [svxyzmat,svid] = gensvalm(refxyz,t);    %Use the almanac-specific
    %                                        %function for sv positions
    prvec = genrng(1,refxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmatref);
    prc=9999*ones(1,24);
    for k = 1:max(size(svid)),
        true_range = norm([svxyzmat(k,:) - refxyz]);
        prc(svid(k)) = prvec(k) - true_range;
    end
%%%User
    clear svid svxyzmat prvec adrvec
    [svxyzmat,svid] = gensvalm(usrxyz,t);
    prvec = genrng(2,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmatusr);
    j=0;
    for k = 1:max(size(svid)),
        if prc(svid(k)) ~= 9999,
           j = j + 1;
           prvec_cr(j) = prvec(k) - prc(svid(k));
           svmat_cr(j,:) = svxyzmat(k,:);
        end
    end
    estusr = olspos(prvec_cr,svmat_cr);
    enuerr(i,:) = ( xyz2enu(estusr(1:3),usrxyz) )';
    terr(i) = estusr(4);  % true clk bias is zero
    clear prvec_cr svmat_cr j k
    waitbar(i/180)
end
close(bar1); 
plot(enuerr(:,1),enuerr(:,2),'*')
axis('square')
axis('equal')
axis([-10 10 -10 10])
grid
title('DGPS {Range Corrections} Positioning Error  -  Yuma-Format Almanac Data')
ylabel('north error (m)')
xlabel('east error (m)')

