%  sndemo22.m      DGPS (Range corrections)
%                  Carrier-Smoothing implemented via Hatch filter
clear all
close all
%    
mpmatusr=mpgen(24,3600,1,22222);
mpmatref=mpgen(24,3600,1,33333);
refllh = [0*pi/180 0*pi/180 0];
refxyz = llh2xyz(refllh);
usrllh = [0.1*pi/180 0.1*pi/180 0];
usrxyz = llh2xyz(usrllh);
loadgps
prrefmat=[]; adrrefmat=[]; prusrmat=[]; adrusrmat=[]; Pvec = [];
i=0;
bar1 = waitbar(0,'Generating DGPS Position...  ');
randn('state',74347098);
for t = 1000:10:2800,
    i=i+1;
%%%Reference Station
    clear svxyzref svidref prref adrref
    [svxyzref,svidref] = gensv(refxyz,t);
    [prref,adrref] = genrng(1,refxyz,svxyzref,svidref,t,...
                           [1 1 0 1 1],[],mpmatref);
    [prsmref,prrefmat,adrrefmat]=...
                   hatch(prref,adrref,svidref,100,prrefmat,adrrefmat);
    prc=9999*ones(1,24);
    for k = 1:max(size(svidref)),
        true_range = norm([svxyzref(k,:) - refxyz]);
        prc(svidref(k)) = prsmref(k) - true_range;
    end
%%%User
    clear svidusr svxyzusr prusr adrusr
    [svxyzusr,svidusr] = gensv(usrxyz,t);
    [prusr,adrusr] = genrng(2,usrxyz,svxyzusr,svidusr,t,...
                            [1 1 0 1 1],[],mpmatusr);
    [prsmusr,prusrmat,adrusrmat]=...
                   hatch(prusr,adrusr,svidusr,100,prusrmat,adrusrmat);
    j=0;
    for k = 1:max(size(svidusr)),
        if prc(svidusr(k)) ~= 9999,
           j = j + 1;
           prvec_cr(j) = prsmusr(k) - prc(svidusr(k));
           svmat_cr(j,:) = svxyzusr(k,:);
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
title('DGPS {Range Corrections} Positioning Error  -  Carrier-Smoothing')
ylabel('north error (m)')
xlabel('east error (m)')

