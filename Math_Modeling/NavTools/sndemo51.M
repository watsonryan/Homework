%  sndemo51.m     Point Positioning with RINEX2 ephemeris and observation files
%                 P-Code pseudoranges are used; Dual-Frequency IONO correction and
%                 standard TROPO correction applied
%
clear all
close all
%
global SQRTSMA
global SVID_MAT TOWSEC PHASE1 PHASE2 C1 P1 P2 D1 D2
global PHASE1LLI PHASE1SS PHASE2LLI PHASE2SS
global C1LLI C1SS P1LLI P1SS P2LLI P2SS
global MARKER_XYZ ANTDELTA OBSINT TIMESTART TIMESTOP CLOCKOFFSET
%
load stkr2581nav
load stkr2581ob

stkrxyz = [678454.119 -4893824.462 4020518.314];    % CORS reference station in Athens, Ohio
usrxyz = stkrxyz;

usr_clk = 0;
sv_clk = zeros(1,32);
c = 299792458;
maskangle = 5;

[r,num_epochs] = size(C1);
bar1 = waitbar(0,'Calculating Position...   ');
n = 0;
for i = 1:5:num_epochs,  %Only processing every 5th sample to speed up this demo
    n = n + 1;   % Counter for the number of samples being processed
    time = TOWSEC(i);     % Time of reception given in GPS time-of-week in seconds
    clear id prvec svxyzmat
    id = find(SVID_MAT(:,i)==1);
    k = 0;
    for j = 1:length(id),  % Loop over all satellites being tracked during this epoch
        if SQRTSMA(id(j)) < 1,
            error('Ephemeris not available for all satellites being tracked')
        end
        if (C1(id(j),i)>SQRTSMA(id(j)))&(P1(id(j),i)>SQRTSMA(id(j)))&(P2(id(j),i)>SQRTSMA(id(j))),     % skip over obviously bad measurements
            tot_est = time -(1/c)*(C1(id(j),i)-usr_clk+sv_clk(id(j)));   % Estimate the time-of-transmission
            [svxyz,E] = svposeph(id(j),tot_est);      % Calculate the position of the satellite
            [sv_clk(id(j)),grpdel] = svclkcorr(id(j),tot_est,E);    % Calculate the satellite clock correction
            if i == 1, estusr(1:3) = usrxyz; end
            svenu = xyz2enu(svxyz,estusr(1:3));     % Convert the satellite position to east-north-up (i.e., local level) coordinates
            el = (180/pi)*atan(svenu(3)/norm(svenu(1:2)));
            if el >= maskangle,
               k = k + 1;       % counter for the satellites being utilized in this epoch
               pL1 = P1(id(j),i) + sv_clk(id(j));    % P-code on L1 pseudorange
               pL2 = P2(id(j),i) + sv_clk(id(j));    % P-code on L2 pseudorange
               prvec(k) = pL1;
               if i > 1,
                  dualcorr = ( (1227.6^2)/(1575.42^2 - 1227.6^2) )*( pL1 - pL2 );  % The dual-freq iono correction
                  tropd = tropocorr(svxyz,estusr(1:3));    % Calculate the tropospheric correction
                  prvec(k) = prvec(k) + dualcorr - tropd;     % Adjust the pseudorange for the iono and tropo corrections
               end
               svxyzr = erotcorr(svxyz,prvec(k));   % Adjust satellite position coordinates for earth rotation correction
               svxyzmat(k,:) = svxyzr';
               svvis(id(j),i) = 1;
            end
        end
    end
    %
    if length(prvec) < 4,
        enuerr(n,:) = NaN;
    else
        estusr = olspos(prvec,svxyzmat);    % Ordinary least-squares position solution
        enuerr(n,:) = ( xyz2enu(estusr(1:3),usrxyz) )';
        usr_clk = estusr(4);
        usr_clk_vec(n) = usr_clk;
    end
    timeofsamples(n) = time;
    %
    waitbar(i/num_epochs,bar1)
end
close(bar1)
close all
plot(enuerr(:,1),enuerr(:,2),'*')
axis('square')
axis('equal')
axis([-10 10 -10 10])
grid
title('GPS P(Y)-Code Point Positioning Error')
ylabel('north error (m)')
xlabel('east error (m)')
text(-7,-5,'TROPO and Dual-Frequency IONO corrections applied')