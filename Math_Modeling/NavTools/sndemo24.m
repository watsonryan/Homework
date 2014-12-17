%  sndemo24.m    Steady-State Differential Carrier-Phase
%                Note effect of atmosphere with lower mask angle
clear all
close all
%
mpmatusr=mpgen(24,3600,1,22222);
mpmatref=mpgen(24,3600,1,33333);
refllh = [0*pi/180 0*pi/180 0];
refxyz = llh2xyz(refllh);
usrllh = [0.1*pi/180 0.1*pi/180 0];
usrxyz = llh2xyz(usrllh);
initref=[0 0 0 0];
initusr=[0 0 0 0];
loadgps
i=0;
bar1 = waitbar(0,'Generating DGPS Position...  ');
randn('state',74347098);
for t = 1000:10:2800,
    i=i+1;
%%%Reference Station
    clear svmatref svmatusr svidref svidusr prref prusr
    [svmatref,svidref] = gensv(refxyz,t,5);
    [trprref,tradrref] = genrng(1,refxyz,svmatref,svidref,t,0); % Error free
    [prref,adrref] = genrng(1,refxyz,svmatref,svidref,t,...
                            [1 1 0 1 1],[],mpmatref);
    estref = olspos(prref,svmatref,initref);  initref=estref;
%%%User
    [svmatusr,svidusr] = gensv(usrxyz,t,5);
    [trprusr,tradrusr] = genrng(2,usrxyz,svmatusr,svidusr,t,0);  % Error free
    [prusr,adrusr] = genrng(2,usrxyz,svmatusr,svidusr,t,...
                           [1 1 0 1 1],[],mpmatusr);
    estusr = olspos(prusr,svmatusr,initusr);  initusr=estusr;
    baspos = mean([estref(1:3); estusr(1:3)]);
    j=0; maxel=0;
    for k = 1:max(size(svidref)),
         id=svidref(k);
         if find(svidusr==id),
            j = j + 1;
            compr(j,1) = prref(k);
            compr(j,2) = prusr(find(svidusr==id));
            compr(j,3) = trprref(k);
            compr(j,4) = trprusr(find(svidusr==id));
            comadr(j,1) = adrref(k);
            comadr(j,2) = adrusr(find(svidusr==id));
            comadr(j,3) = tradrref(k);
            comadr(j,4) = tradrusr(find(svidusr==id));
            comsv(j,:) = svmatref(k,:);
            vec = comsv(j,:) - baspos; 
            unit_vec(j,:) = vec / norm(vec);
            enu = xyz2enu(comsv(j,:),baspos);
            el(j) = atan2(enu(3),norm(enu(1:2)));
            if el(j) > maxel,
               maxel = el(j);
               keyid = j;
            end
         end
    end
    numcomsv = j;
%
    clear sd_pr dd_pr H dd_adr trsd_pr sd_adr trsd_adr tr_N_lam
    clear trdd_pr trdd_adr
    sd_pr = compr(:,2) - compr(:,1);
    trsd_pr = compr(:,4) - compr(:,3);
    sd_adr = comadr(:,2) - comadr(:,1);
    trsd_adr = comadr(:,4) - comadr(:,3);
    j = 0;
    for k = 1:numcomsv,
        if k ~= keyid,
           j = j + 1;
           dd_pr(j,1) = sd_pr(keyid) - sd_pr(k);
           trdd_pr(j,1) = trsd_pr(keyid) - trsd_pr(k);
           dd_adr(j,1) = sd_adr(keyid) - sd_adr(k);
           trdd_adr(j,1) = trsd_adr(keyid) - trsd_adr(k);
           H(j,:) = unit_vec(keyid,:) - unit_vec(k,:);
        end
    end
%
    base_est = H\dd_pr;
    ddestusr = refxyz - base_est';
    err_pr(i,:) = ddestusr - usrxyz;
%
    tr_N_lam = trdd_adr - trdd_pr;
%
    base_est = H\(dd_adr - tr_N_lam);
    ddestusr = refxyz - base_est';
%%%    err_adr(i,:) = ddestusr - usrxyz;
    enuerr(i,:) = ( xyz2enu(ddestusr,usrxyz) )';
%
    clear sd_pr dd_pr H dd_adr trsd_pr sd_adr trsd_adr tr_N_lambda
    clear trdd_pr trdd_adr ddestusr base_est
    clear svidref svidusr svmatref svmatusr t tradrref tradrusr
    clear trprref trprusr unit_vec vec U adrref adrusr c comadr compr 
    clear comsv el enu estref estusr id j k keyid maxel numcomsv 
    clear prref prusr r
    waitbar(i/180)
end
close(bar1); 
plot(enuerr(:,1),enuerr(:,2),'*')
axis('square')
axis('equal')
axis([-.5 .5 -.5 .5])
grid
title('Steady-State Differential Carrier-Phase Positioning Error')
ylabel('north error (m)')
xlabel('east error (m)')
