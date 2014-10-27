%
%   insgps_18_with_jump.m
%
%   Loosely-Coupled 18-State INS/GPS Kalman Filter
%
%   This program is the same as INSGPS_18_STATE but
%   here we purposely reduce the number of visible
%   satellites (towards the end of the run) to show
%   the filters response to a change in the GPS
%   constellation.
%
%   The state vector is:
%    x1 = delta x position
%    x2 = delta y position
%    x3 = delta z position
%    x4 = delta x velocity
%    x5 = delta y velocity
%    x6 = delta z velocity
%    x7 = psi-angle x component
%    x8 = psi-angle y component
%    x9 = psi-angle z component
%   x10 = body-x accel bias
%   x11 = body-y accel bias
%   x12 = body-z accel bias
%   x13 = body-x gyro bias
%   x14 = body-y gyro bias
%   x15 = body-z gyro bias
%   x16 = x-component of bias in GPS estimated position
%   x17 = y-component of bias in GPS estimated position
%   x18 = z-component of bias in GPS estimated position
%
%    Note:  the position and velocity errors are expressed in the
%           so-called local-level frame or simply L (or l) frame.
%           In the L-frame:  x=east; y=north; z=up
%
%  September 2005
%  Copyright (c) 2005 by GPSoft LLC
%  All Rights Reserved.
%
%  Revision History
%  Nov 2008:  Fixed bug with gyro bias states in F matrix

clear all
close all

randn('state',0)

mpmat=mpgen(24,3600,1,54321);  % simulate GPS multipath error
loadgps  % load the GPS constellation parameters into global memory
gpstime = 40000;  % specify the GPS time-of-day for the simulation

load dyn_flt_ins_dat    % output from the program DYN_FLIGHT_INS.M

C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

H = zeros(3,18);  H(1,1) = 1;  H(2,2) = 1;  H(3,3) = 1;
%                                 % The measurement vector z consists of
%                                 % an x, y and z position differences only.
%                                 % Since we are working in the L frame
%                                 % this corresponds to east, north and up
%                                 % position differences (since alpha=0).

%  Defining the components of the initial prediction error covariance
%  matrix
init_pos_rms = 3;  % initial position error (RMS in meters)
Peast_pos = init_pos_rms^2;  Pnorth_pos = init_pos_rms^2;  Pup_pos = 0;
Peast_vel = 2^2;  Pnorth_vel = 2^2;  Pup_vel = 0;
Ppsi_x = 0.0001^2;  Ppsi_y = 0.0001^2;  Ppsi_z = 0.0001^2;
Pacc_x = (100*9.81e-6)^2; Pacc_y = (100*9.81e-6)^2; Pacc_z = (100*9.81e-6)^2;
Pgyr_x = (0.05)^2; Pgyr_y = (0.05)^2; Pgyr_z = (0.05)^2;

P_pre = zeros(18,18);
P_pre(1,1)=Peast_pos; P_pre(2,2)=Pnorth_pos; P_pre(3,3)=Pup_pos;
P_pre(4,4)=Peast_vel; P_pre(5,5)=Pnorth_vel; P_pre(6,6)=Pup_vel;
P_pre(7,7)=Ppsi_x; P_pre(8,8)=Ppsi_y; P_pre(9,9)=Ppsi_z;
P_pre(10,10)=Pacc_x; P_pre(11,11)=Pacc_y; P_pre(12,12)=Pacc_z;
P_pre(13,13)=Pgyr_x; P_pre(14,14)=Pgyr_y; P_pre(15,15)=Pgyr_z;
P_pre(16,16)=init_pos_rms^2; P_pre(17,17)=init_pos_rms^2; P_pre(18,18)=init_pos_rms^2;

P_pre_KF = P_pre;   % initial prediction error covariance matrix

%  measurement error covariance matrix
R = init_pos_rms^2*eye(3);

%  In this section we are defining the matrices which are used to 
%  numerically calculate the system noise covariance matrix (Q).
%  The two 'sigma' numbers are rough estimates of the standard deviation of
%  the white noise input to the first-order Gauss-Markov processes which
%  are being used to model the accelerometer and gyro biases
sigma_acc = 0.3*9.81e-6;
sigma_gyr = 1e-9;
G = zeros(18,18); 
G(10,10)=sigma_acc; G(11,11)=sigma_acc; G(12,12)=sigma_acc; 
G(13,13)=sigma_gyr; G(14,14)=sigma_gyr; G(15,15)=sigma_gyr; 
W = zeros(18,18); 
W(10,10)=1; W(11,11)=1; W(12,12)=1;
W(13,13)=1; W(14,14)=1; W(15,15)=1;

%  Time constants (in seconds) for the first-order Gauss-Markov processes
%  being used to model the accelerometer and gyro biases
tau_accel = 100000;
tau_gyro = 100000;
update = 0;  count = 0; i_KF = 0;

X_pre = zeros(18,1);   % Initialize the state prediction

est_lat_KF(1) = lat_prof(1);   % INITIALIZATION IN THIS SECTION
est_lon_KF(1) = lon_prof(1);
est_alpha_KF(1) = 0;
est_vel_l_KF(1,:) =  vel_l(1,:);
est_roll_KF(1) = est_roll(1);
est_pitch_KF(1) = est_pitch(1);
est_yaw_KF(1) = est_yaw(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,' Time Loop ');
for i = 2:npts,

    % Forming the system dynamics matrix
    [rm,rp] = radicurv(est_lat(i));  radius_e = sqrt(rm*rp);
    accel_vect_L = del_Vl(:,i)*(1/tdint(i));
    F11 = -1*skewsymm(omega_el_L(:,i));    % L-frame
    F12 = eye(3);
    F13 = 0*eye(3);
    F21=eye(3); F21(1,1)=-g_extr(i)/radius_e; 
    F21(2,2)=-g_extr(i)/radius_e; F21(3,3)=2*g_extr(i)/(radius_e+est_height(i,1));
    F22 = -1*skewsymm(2*omega_ie_L(:,i) + omega_el_L(:,i));   % L-frame
    F23 = skewsymm(accel_vect_L);
    F31 = 0*eye(3);
    F32 = 0*eye(3);
    F33 = (-1)*skewsymm(omega_ie_L(:,i) + omega_el_L(:,i));   % L-frame
    F_ins = [F11 F12 F13; F21 F22 F23; F31 F32 F33];

    F = zeros(18,18); F(1:9,1:9) = F_ins; 
    F(10:12,10:12) = (-1/tau_accel)*eye(3);  F(4:6,10:12) = C*est_DCMbn(:,:,i);
    F(13:15,13:15) = (-1/tau_gyro)*eye(3);  F(7:9,13:15) = -1*C*est_DCMbn(:,:,i);

    
    %  Numerical calculation of the state transition matrix (PHI) and
    %  system noise covariance matrix (Q)
    A = zeros(36,36);
    A(1:18,1:18) = -1*F;
    A(1:18,19:36) = G*W*G';
    A(19:36,19:36) = F';
    A = A*tdint(i);

    B = expm(A);

    PHI_trans = B(19:36,19:36);
    PHI = PHI_trans';

    Q = PHI*B(1:18,19:36);

    % Calculate the Kalman gain
    K = P_pre_KF*H'*inv(H*P_pre_KF*H' + R);
    
    % keep track of elapsed time (in seconds) since last Kalman update
   count = count + tdint(i);

    %%if count >= 630,   % 10.5 minute intervals (i.e., 1/8 Schuler cycle)
    %%if count >= 60,  % 1 minute intervals
    if count >= 1,    % 1 Hz updates
        update = 1;
        update_interval = count;
        count = 0;
    end
    
    if update == 1,
        i_KF = i_KF + 1;
        %%  Simulate GPS Receiver output here
        gpstime = gpstime + update_interval;
        tru_pos_xyz = llh2xyz([lat_prof(i) lon_prof(i) height_prof(i)]);
        [svxyzmat,svid] = gensv(tru_pos_xyz,gpstime,5);
        %%[prvec,adrvec] = genrng(1,tru_pos_xyz,svxyzmat,svid,gpstime,[1 1 0 1 1],[],mpmat);
        [prvec,adrvec] = genrng(1,tru_pos_xyz,svxyzmat,svid,gpstime,[1 0 0 0 0],[],mpmat);

        if i_KF > 1020,        % artificially reducing the number of 'visible' satellites
            numsv = 4;
        else
            numsv = length(prvec);
        end
        
        gps_xyzt = olspos(prvec(1:numsv),svxyzmat(1:numsv,:));
        gps_llh = xyz2llh(gps_xyzt(1:3));

        ins_lat = lat_prof(i) + est_lat_err(i);
        ins_lon = lon_prof(i) + est_lon_err(i);

        %  Kalman observation vector (Z) formed by position differences.
        %  Angular position differences converted to linear differences.
        Z(1,1) = (ins_lon - gps_llh(2))*rp*cos(ins_lat);  % east component
        Z(2,1) = (ins_lat - gps_llh(1))*rm;               % north component
        Z(3,1) = height_prof(i) - gps_llh(3);             % up component

        X_est = X_pre + K*(Z - H*X_pre);

        P_est_KF = (eye(18) - K*H)*P_pre_KF;

        update = 0;
    else
        X_est = X_pre;
            X_est(3) = 0;  % Clamp the vertical error during coasting
            X_est(6) = 0;
        P_est_KF = P_pre_KF;
    end
    X_pre = PHI*X_est;        
    P_pre_KF = PHI*P_est_KF*PHI' + Q;
    
    state_est(:,i) = X_est;
    x_rms_KF(i) = sqrt( P_est_KF(1,1) );
    y_rms_KF(i) = sqrt( P_est_KF(2,2) );
    x_v_rms_KF(i) = sqrt( P_est_KF(4,4) );
    y_v_rms_KF(i) = sqrt( P_est_KF(5,5) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta(1,1) = -X_est(2)/radius_e;   % converting position error from linear
    theta(2,1) = X_est(1)/radius_e;    % units to angular units
    theta(3,1) = tan(est_lat(i))*theta(2);  % this is for alpha=0 and ENU frame;

    psi = X_est(7:9);
    phi_angle = psi + theta;         % computing total attitude error

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % This is the Kalman estimate of the true DCM
    est_DCMbn_KF = C*(eye(3) + skewsymm(phi_angle))*C*est_DCMbn(:,:,i); 
    
    eul_vect = dcm2eulr(est_DCMbn_KF);
    est_roll_KF(i) = eul_vect(1);
    est_pitch_KF(i) = eul_vect(2);
    est_yaw_KF(i) = eul_vect(3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % This is the Kalman estimate of the true DCM
    est_DCMel_KF = (eye(3) + skewsymm(theta))*est_DCMel(:,:,i);
    
    llw_vect = dcm2llw(est_DCMel_KF);
    est_lat_KF(i) = llw_vect(1);
    est_lon_KF(i) = llw_vect(2);
    est_alpha_KF(i) = llw_vect(3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    est_vel_l_KF(i,:) = ( (eye(3) + skewsymm(theta))*(vel_l(i,:)' - X_est(4:6)) )';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    waitbar(i/npts,h)
    
    totnumsv(i) = numsv;
end
close(h)


h = waitbar(0,' Computing Position, Velocity, Attitude Errors ');
N = max(size(est_lat0));                           % Compute horizontal position
for i = 1:N,                                       % error by finding the ENU
   truxyz = llh2xyz([est_lat0(i) est_lon0(i) 0]);  % coordinates of the 
   insxyz_unaided = llh2xyz([est_lat(i) est_lon(i) 0]);    % INS-derived position 
   enu = xyz2enu(insxyz_unaided,truxyz);                   % relative to the truth
   east_pos_err_unaided(i) = enu(1);
   north_pos_err_unaided(i) = enu(2);
   up_pos_err_unaided(i) = enu(3);
   horz_pos_err_unaided(i) = norm(enu(1:2));

   insxyz_KF = llh2xyz([est_lat_KF(i) est_lon_KF(i) 0]);
   enu = xyz2enu(insxyz_KF,truxyz);
   east_pos_err_KF(i) = enu(1);
   north_pos_err_KF(i) = enu(2);
   up_pos_err_KF(i) = enu(3);
   horz_pos_err_KF(i) = norm(enu(1:2));

   waitbar(i/N,h)
end
close(h)

lat_err_unaided = est_lat-est_lat0;
lon_err_unaided = est_lon-est_lon0;
alpha_err_unaided = est_alpha - est_alpha0;
vel_l_err_unaided = vel_l - vel_l_0;
roll_err_unaided = est_roll - est_roll0;
pitch_err_unaided = est_pitch - est_pitch0;
yaw_err_unaided = est_yaw - est_yaw0;

lat_err_KF = est_lat_KF-est_lat0;
lon_err_KF = est_lon_KF-est_lon0;
alpha_err_KF = est_alpha_KF - est_alpha0;
vel_l_err_KF = est_vel_l_KF - vel_l_0;
roll_err_KF = est_roll_KF - est_roll0;
pitch_err_KF = est_pitch_KF - est_pitch0;
yaw_err_KF = est_yaw_KF - est_yaw0;


close all

subplot(111)
plot(time/60,lat_err_unaided*180/pi,time/60,lon_err_unaided*180/pi,...
     time/60,lat_err_KF*180/pi,time/60,lon_err_KF*180/pi)
title('Loosely-Coupled 18-State INS/GPS Filter')
ylabel('error in degrees')
xlabel('time in minutes')
legend('lat err unaided', 'long err unaided', 'lat err KF', 'lon err KF',0)

figure
subplot(111)
plot(time/60,vel_l_err_unaided(:,1:2),time/60,vel_l_err_KF(:,1:2))
title('Loosely-Coupled 18-State INS/GPS Filter')
ylabel('velocity error in m/s')
xlabel('time in minutes')
legend('east unaided','north unaided', 'east KF', 'north KF',0)

figure
subplot(111)
plot(time/60,roll_err_unaided*180/pi,time/60,pitch_err_unaided*180/pi,...
     time/60,roll_err_KF*180/pi,time/60,pitch_err_KF*180/pi)
title('Loosely-Coupled 18-State INS/GPS Filter - Euler Angle Errors')
ylabel('error in degrees')
xlabel('time in minutes')
legend('roll unaided','pitch unaided', 'roll KF', 'pitch KF',0)

figure
subplot(111)
plot(time/60,yaw_err_unaided*180/pi,time/60,yaw_err_KF*180/pi)
title('Loosely-Coupled 18-State INS/GPS Filter')
ylabel('yaw error in degrees')
xlabel('time in minutes')
legend('yaw unaided', 'yaw KF',0)

figure
plot(time/60,horz_pos_err_unaided,time/60,horz_pos_err_KF)
title('Loosely-Coupled 18-State INS/GPS Filter - Horizontal Position Error')
ylabel('error in meters')
xlabel('time in minutes')
legend('unaided','KF',0)

figure
plot(time/60,horz_pos_err_KF)
title('Loosely-Coupled 18-State INS/GPS Filter - Horizontal Position Error')
ylabel('error in meters')
xlabel('time in minutes')

figure
plot(time/60,x_rms_KF,time/60,y_rms_KF)
title('Loosely-Coupled 18-State INS/GPS Filter - Covariance Analysis')
ylabel('position rms error in meters')
xlabel('time in minutes')
legend('x','y',0)

figure
plot(time/60,x_v_rms_KF,time/60,y_v_rms_KF)
title('Loosely-Coupled 18-State INS/GPS Filter Covariance Analysis')
ylabel('velocity rms error in meters/sec')
xlabel('time in minutes')
legend('x','y',0)


figure
subplot(211)
plot(time/60,horz_pos_err_KF)
title('Tight-Coupling and Degraded GPS Coverage')
ylabel('horizontal position error in meters')
subplot(212)
plot(time/60,totnumsv)
xlabel('time in minutes')
ylabel('number of satellites')
