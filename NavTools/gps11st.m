%  gps11st.m
%
%  Example of an 11-State GPS Extended Kalman Filter suitable
%  for moderate-dynamic users
%
%  The state vector is:
%   x1 = delta x position
%   x2 = delta x velocity
%   x3 = delta x acceleration
%   x4 = delta y position
%   x5 = delta y velocity
%   x6 = delta y acceleration
%   x7 = delta z position
%   x8 = delta z velocity
%   x9 = delta z acceleration
%   x10 = clock offset
%   x11 = clock drift
%

%	References:
%
%	Brown, R. G. and P. Y. C. Hwang, "Introduction to Random
%	Signals and Applied Kalman Filtering," 3rd edition, John Wiley & Sons,
%	1997.
%
%   Levy, L., "Integration of GPS with Inertial Navigation Systems,"
%   short course notes, Navtech Seminars, Springfield, VA, 2003.
%
%  March 2004; August 2005
%  Copyright (c) 2004-2005 by GPSoft LLC
%  All Rights Reserved.
%
clear all
close all
randn('state',0)
%
load dogleg1hz     % truth data are generated in an
%                  % East-North-Up coordinate system
usrenu = tru_pos_L;
%
mpmat=mpgen(24,3600,1,54321);

% Specify the origin of the local-level coordinate frame
lat = 39; lon = -81;
orgllh = [39*pi/180 -81*pi/180 0];
orgxyz = llh2xyz(orgllh);

% Now compute the direction cosine matrix relating ECEF
% xyz coordinates to a local Up-East-North frame
clat = cos(-lat*pi/180); clon = cos(lon*pi/180);
slat = sin(-lat*pi/180); slon = sin(lon*pi/180);
Cz = [clon slon 0; -slon clon 0; 0 0 1];
Cy = [clat 0 -slat; 0 1 0; slat 0 clat];
Cxyz2uen=Cy*Cz;

loadgps
speedoflight = 299792458;

deltat=1;   % specify Kalman update interval

% generate system noise covariance matrix and state transition matrix
Sp=10; Sf=1e-3; Sg=1e-3; tauaccel=60;
[Q,PHI] = q_gen_11_gps(deltat,Sp,Sf,Sg,tauaccel);

sigmaz = 1;  % standard deviation of pseudorange measurement
%            % noise in meters
I = eye(11);
x_pre = zeros(11,1);   % initial state vector prediction
P_pre = 10*eye(11);    % initial prediction error covariance matrix

EndLoop = max(size(usrenu));
bar1 = waitbar(0,'Generating User Solutions...  ');
for i = 1:EndLoop,
    t = tru_time(i);
    usrxyz=enu2xyz(usrenu(i,:),orgxyz);   % convert current true user position
    % from east-north-up to ECEF x-y-z coordinates
    
    [svxyzmat,svid] = gensv(usrxyz,t,5);  % determine current satellites in view
    
    %[prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 1 0 1 1],[],mpmat);
    [prvec,adrvec] = genrng(1,usrxyz,svxyzmat,svid,t,[1 0 0 0 0],[],mpmat);

    % simulate the receiver clock drift as 100 ns/s and an initial offset
    % of 10 milliseconds
    rx_clk_offset(i) = 10e-3*speedoflight + speedoflight*i*1e-7;
    
    prvec = prvec + rx_clk_offset(i);
    adrvec = adrvec + rx_clk_offset(i);
    
    if i == 1,    % if this is the first data point, simply perform ordinary least squares positioning
        estusr = olspos(prvec,svxyzmat);
        x_nom=zeros(11,1); x_nom(1)=estusr(1); x_nom(4)=estusr(2);
                          x_nom(7)=estusr(3); x_nom(10)=estusr(4);
    else
        n = length(svid);
        nom_pos = [x_nom(1) x_nom(4) x_nom(7)];
        Hpartial = hmat(svxyzmat,nom_pos);        % calculate the direction cosine matrix elements
        H = zeros(n,11); H(:,1)=Hpartial(:,1); H(:,4)=Hpartial(:,2);
        H(:,7)=Hpartial(:,3); H(:,10)=Hpartial(:,4);
        for j = 1:n,
            pr0 = norm(svxyzmat(j,:) - nom_pos) + x_nom(10);  % form the calculated pseudorange to the nominal position
            z(j,1) = prvec(j) - pr0;        % form the Kalman observation (measurement)
        end

        R = sigmaz*eye(length(prvec));   % form the measurement error covariance matrix
        
        K = P_pre*H'*inv(H*P_pre*H' + R);   % form the Kalman gain matrix
        
        x_est = x_pre + K*(z - H*x_pre);  % form the Kalman estimate
        
        P = (I - K*H)*P_pre;   % form the Kalman estimation error covariance matrix
        
        x_tot_est = x_nom + x_est;   % since the estimated state vector is composed of delta-quantities,
        %  the total estimate is formed by summing the nominal with the
        %  estimated delta quantities
             estusr = [x_tot_est(1) x_tot_est(4) x_tot_est(7) x_tot_est(10)];
             estusr_rate = [x_tot_est(2) x_tot_est(5) x_tot_est(8) x_tot_est(11)];
             estuenvel=Cxyz2uen*estusr_rate(1:3)';
             estenuvel(i,:) = [estuenvel(2) estuenvel(3) estuenvel(1)];
        
        x_nom = PHI*x_tot_est;   % since this is an extended Kalman filter, we use the current
        % estimate to form the next nominal point (i.e., the next point
        % around which the linearization will take place
        x_pre = zeros(11,1);  % again, since this is an EKF, the state prediction is zeros
        P_pre = PHI*P*PHI' + Q;  % form the prediction error covariance matrix
    end    % end 'if i==1' loop

    estenu(i,:) = ( xyz2enu(estusr(1:3),orgxyz) )';  % convert position from ECEF x-y-z to ENU
    estclockb(i) = estusr(4);
    err(i,1:3) = estenu(i,1:3) - usrenu(i,:);  % calculate position error in east-north-up
    terr(i) = estusr(4) - rx_clk_offset(i);  % calculate timing error (i.e., receiver clock bias estimation error)
    waitbar(i/EndLoop)   
end    % end 'for i=1:EndLoop' loop
close(bar1);

plot(usrenu(:,1),usrenu(:,2),'-',estenu(:,1),estenu(:,2),'*')
title('True and Estimated Trajectories')
ylabel('north direction [meters]')
xlabel('east direction [meters]')

figure
plot(tru_time,err(:,1),tru_time,err(:,2),tru_time,err(:,3),tru_time,terr)
ylabel('position error in meters')
xlabel('run time in seconds')
title('GPS 11-State Kalman Estimation Errors')
legend('x','y','z','b',0)

figure
plot(tru_time,err(:,1),tru_time,err(:,2),tru_time,err(:,3),tru_time,terr)
ylabel('position error in meters')
xlabel('run time in seconds')
title('ZOOM VIEW: GPS 11-State EKF Estimation Errors')
legend('x','y','z','b',0)
axis([0 150 -7 3])

figure
vel_err = estenuvel - tru_vel_L;
plot(tru_time,vel_err(:,1),tru_time,vel_err(:,2),tru_time,vel_err(:,3))
ylabel('velocity error in meters/second')
xlabel('run time in seconds')
title('GPS 11-State Kalman Estimation Errors')
legend('x','y','z',0)

figure
vel_err = estenuvel - tru_vel_L;
plot(tru_time,vel_err(:,1),tru_time,vel_err(:,2),tru_time,vel_err(:,3))
ylabel('velocity error in meters/second')
xlabel('run time in seconds')
title('ZOOM VIEW: GPS 11-State EKF Estimation Errors')
legend('x','y','z',0)
axis([0 150 -6 6])
