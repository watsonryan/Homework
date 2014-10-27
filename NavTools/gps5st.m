%  gps5st.m
%
%  Example of a 5-State GPS Extended Kalman Filter suitable
%  for stationary users
%
%  The state vector is:
%   x1 = delta x position
%   x2 = delta y position
%   x3 = delta z position
%   x4 = delta clock offset
%   x5 = delta clock drift
%

%	References:
%   Brown, R. G. and P. Y. C. Hwang, "Introduction to Random Signals and
%   Applied Kalman Filtering," 3rd edition, John Wiley & Sons, New York,
%   1997.
%
%   Levy, L., "Integration of GPS with Inertial Navigation Systems,"
%             short course notes, Navtech Seminars, Springfield, VA, 2003.
%
%  August 2005
%  Copyright (c) 2005 by GPSoft LLC
%  All Rights Reserved.
%
clear all
close all
randn('state',0)
%
mpmat=mpgen(24,3600,1,54321);
%
% Specify the true user position
lat = 39; lon = -81;        
tru_usrllh = [39*pi/180 -81*pi/180 0];
tru_usrxyz = llh2xyz(tru_usrllh);

loadgps  % load the satellite orbital parameters into Global memory
speedoflight = 299792458;

deltat=1;        % specify Kalman update interval
% generate system noise covariance matrix and state transition matrix
Sp = 0.1; Sf = 0.1; Sg = 0.1;
[Q,PHI] = q_gen_5_gps(deltat,Sp,Sf,Sg);

sigmaz = 1;  % standard deviation of pseudorange measurement
             % noise in meters

I = eye(5);

x_pre = zeros(5,1);  % initial state prediction
P_pre = 50*eye(5);   % initial prediction error covariance matrix
%                    % this is quite conservative since the initial
%                    % point of linearization is formed using the ordinary
%                    % least squares (OLS) solution

Ntot = 1000;  % length of simulation
bar1 = waitbar(0,'Generating User Solutions...  ');
for i = 1:Ntot,
    time(i) = 42000 + i;
    
    [svxyzmat,svid] = gensv(tru_usrxyz,time(i),5);
    
    %   The following simulates full thermal noise, tropo delay, multipath
    %   and iono delay (i.e., no corrections)
    %[prvec,adrvec] = genrng(1,tru_usrxyz,svxyzmat,svid,time(i),[1 1 0 1 1],[],mpmat);
    
    %   The following simulates an ideal differential-correction scenario
    %   where all errors (except noise) have been eliminated
    %[prvec,adrvec] = genrng(1,tru_usrxyz,svxyzmat,svid,time(i),[1 0 0 0 0],[],mpmat);
    
    %   The following simulates a civilian single-frequency user employing
    %   the broadcast iono correction and standard tropo correction models.
    %   The standard tropo correction removes 90% of the raw tropo delay
    %   and the broadcast iono correction removes, on average, 50% of the
    %   raw iono delay
    [prvec,adrvec] = genrng(1,tru_usrxyz,svxyzmat,svid,time(i),[1 0.1 0 1 0.5],[],mpmat);

    % simulate the receiver clock drift as 100 ns/s and an initial offset
    % of 10 milliseconds
    rx_clk_offset(i) = 10e-3*speedoflight + speedoflight*i*1e-7;
        prvec = prvec + rx_clk_offset(i);
        adrvec = adrvec + rx_clk_offset(i);

    if i == 1,    % if this is the first data point, simply perform ordinary least squares positioning
        estusr = olspos(prvec,svxyzmat);
        x_nom=zeros(5,1); x_nom(1)=estusr(1); x_nom(2)=estusr(2);
                          x_nom(3)=estusr(3); x_nom(4)=estusr(4);
    else
        n = length(svid);
        nom_pos = [x_nom(1) x_nom(2) x_nom(3)];
        Hpartial = hmat(svxyzmat,nom_pos);        % calculate the direction cosine elements
        H = zeros(n,5); H(:,1:4)=Hpartial;
        for j = 1:n,
            pr0 = norm(svxyzmat(j,:) - nom_pos) + x_nom(4);  % form the calculated pseudorange to the nominal position
            z(j,1) = prvec(j) - pr0;        % form the Kalman observation (measurement)
        end

        R = sigmaz*eye(length(prvec));   % form the measurement error covariance matrix

        K = P_pre*H'*inv(H*P_pre*H' + R);   % form the Kalman gain matrix

        x_est = x_pre + K*(z - H*x_pre);  % form the Kalman estimate
        P = (I - K*H)*P_pre;   % form the Kalman estimation error covariance matrix

        x_tot_est = x_nom + x_est;   % since the estimated state vector is composed of delta-quantities,
        %  the total estimate is formed by summing the nominal with the
        %  estimated delta quantities
        estusr = [x_tot_est(1) x_tot_est(2) x_tot_est(3) x_tot_est(4)];
        
        x_nom = PHI*x_tot_est;   % since this is an extended Kalman filter, we use the current
        % estimate to form the next nominal point (i.e., the next point
        % around which the linearization will take place
        x_pre = zeros(5,1);  % again, since this is an EKF, the state prediction is zeros
        P_pre = PHI*P*PHI' + Q;  % form the prediction error covariance matrix

    end    % end 'if i==1' loop
    
    estenu_err(i,:) = ( xyz2enu(estusr(1:3),tru_usrxyz) )';  % calculate position error in east-north-up
    estclockb(i) = estusr(4);
    terr(i) = estusr(4) - rx_clk_offset(i);  % calculate timing error (i.e., receiver clock bias estimation error)
    
    ols_est = olspos(prvec,svxyzmat);   % same as above but this is for the OLS solution
    ols_enu_err(i,:) = ( xyz2enu(ols_est(1:3),tru_usrxyz) )';
    ols_terr(i) = ols_est(4) - rx_clk_offset(i);

    waitbar(i/Ntot)   

end   % end 'for i=1:Ntot' loop
close(bar1);

figure
plot(time,estenu_err(:,1),time,estenu_err(:,2),time,estenu_err(:,3),time,terr)
axis([42000 43000 -4 16])
ylabel('position error in meters')
xlabel('GPS time of day in seconds')
title('GPS 5-State Kalman Estimation Errors')
legend('x','y','z','b',0)

figure
plot(time,ols_enu_err(:,1),time,ols_enu_err(:,2),time,ols_enu_err(:,3),time,ols_terr)
axis([42000 43000 -4 16])
ylabel('position error in meters')
xlabel('GPS time of day in seconds')
title('GPS OLS Estimation Errors')
legend('x','y','z','b',0)

