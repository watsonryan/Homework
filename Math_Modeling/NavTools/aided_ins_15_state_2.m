%
%   aided_ins_15_state_2.m
%
%   Similar to AIDED_INS_15_STATE.M but here we are using a
%   different and longer flight path (see LOAD command below).
%   In addition, the Q matrix elements are increased to keep
%   the filter stable (i.e., if we used the same Q matrix as
%   the one in AIDED_INS_15_STATE.M then the filter would diverge).
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

fly_over_rms = 3;  % RMS error (in meters) of the external 2-D position
%                  % fix (in each dimension); also known as a fly-over fix

load dyn_flt_2_ins_dat    % output from the program DYN_FLIGHT_2_INS.M

randn('state',0)

C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

H = zeros(2,15);  H(1,1) = 1;  H(2,2) = 1;
%                                 % The measurement vector z consists of
%                                 % an x and y position difference only.
%                                 % Since we are working in the L frame
%                                 % this corresponds to east and north
%                                 % position difference (since alpha=0).

Peast_pos = fly_over_rms^2;  Pnorth_pos = fly_over_rms^2;  Pup_pos = 0;
Peast_vel = 2^2;  Pnorth_vel = 2^2;  Pup_vel = 0;
Ppsi_x = 0.0001^2;  Ppsi_y = 0.0001^2;  Ppsi_z = 0.0001^2;
Pacc_x = (100*9.81e-6)^2; Pacc_y = (100*9.81e-6)^2; Pacc_z = (100*9.81e-6)^2;
Pgyr_x = (0.05)^2; Pgyr_y = (0.05)^2; Pgyr_z = (0.05)^2;

P_pre = zeros(15,15);
P_pre(1,1)=Peast_pos; P_pre(2,2)=Pnorth_pos; P_pre(3,3)=Pup_pos;
P_pre(4,4)=Peast_vel; P_pre(5,5)=Pnorth_vel; P_pre(6,6)=Pup_vel;
P_pre(7,7)=Ppsi_x; P_pre(8,8)=Ppsi_y; P_pre(9,9)=Ppsi_z;
P_pre(10,10)=Pacc_x; P_pre(11,11)=Pacc_y; P_pre(12,12)=Pacc_z;
P_pre(13,13)=Pgyr_x; P_pre(14,14)=Pgyr_y; P_pre(15,15)=Pgyr_z;

P_pre_KF = P_pre;   % initial prediction error covariance matrix

R = [fly_over_rms^2 0; 0 fly_over_rms^2];  % measurement error covariance matrix

sigma_acc = 0.3*9.81e-6;
sigma_gyr = 1e-9;
G = zeros(15,15); 
G(10,10)=sigma_acc; G(11,11)=sigma_acc; G(12,12)=sigma_acc; 
G(13,13)=sigma_gyr; G(14,14)=sigma_gyr; G(15,15)=sigma_gyr; 
W = zeros(15,15); 
W(10,10)=1; W(11,11)=1; W(12,12)=1;
W(13,13)=1; W(14,14)=1; W(15,15)=1;
%%%%%%%%%%%%%%%%%%%%%%% NOTE non-zero values for attitude states added here
        G(7:9,7:9) = 1e-4*eye(3);
        W(7:9,7:9) = eye(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_accel = 100000;
tau_gyro = 100000;
update = 0;  count = 0;

X_pre = zeros(15,1);

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
    
    % Generate the inertial error model system dynamics matrix
    F_ins = F_ins_gen(est_lat(i),del_Vl(:,i),tdint(i),...
        omega_el_L(:,i),omega_ie_L(:,i),g_extr(i),est_height(i,1));

    % Generate the system dynamics matrix for the whole filter
    F = zeros(15,15); F(1:9,1:9) = F_ins; 
    F(10:12,10:12) = (-1/tau_accel)*eye(3);  F(4:6,10:12) = C*est_DCMbn(:,:,i);
    F(13:15,13:15) = (-1/tau_gyro)*eye(3);  F(7:9,13:15) = -1*C*est_DCMbn(:,:,i);
    
    % Generate the system noise covariance matrix and the state transition
    % matrix
    [Q,PHI] = q_gen(tdint(i),F,G,W,15);

    % Calculate the Kalman gain
    K = P_pre_KF*H'*inv(H*P_pre_KF*H' + R);
    
   % keep track of elapsed time (in seconds) since last Kalman update
   count = count + tdint(i);

    %%if count >= 630,   % 10.5 minute intervals (i.e., 1/8 Schuler cycle)
    %%if count >= 60,  % 1 minute intervals
    if count >= 1,    % 1 Hz updates
        update = 1;
        count = 0;
    end
    
    if update == 1,
        % Generate the free-inertial estimates of lat, long and ECEF x-y-z
        ins_lat = lat_prof(i) + est_lat_err(i);
        ins_lon = lon_prof(i) + est_lon_err(i);
        est_pos_xyz = llh2xyz([ins_lat ins_lon height_prof(i)]);

        % Calcuate the ECEF x-y-z coordinates of the true position
        tru_pos_xyz = llh2xyz([lat_prof(i) lon_prof(i) height_prof(i)]);

        % Calculate the free-inertial position error in east-north-up
        % coordinates
        ins_pos_err_enu = xyz2enu(est_pos_xyz,tru_pos_xyz);

        % Simulate the errors of the position aiding source
        fly_over_fix_errors = fly_over_rms*randn(2,1);

        % Simulate the position difference between INS and position aiding
        % source
        Z = ins_pos_err_enu(1:2,1) + fly_over_fix_errors;  
        
        % Compute Kalman estimate
        X_est = X_pre + K*(Z - H*X_pre);
            X_est(3) = 0;
            X_est(6) = 0;    % Assuming perfect vertical aiding
        P_est_KF = (eye(15) - K*H)*P_pre_KF;

        update = 0;
    else
        X_est = X_pre;
            X_est(3) = 0;
            X_est(6) = 0;    % Assuming perfect vertical aiding
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
    [rm,rp] = radicurv(est_lat(i));  radius_e = sqrt(rm*rp);
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
title('15-State Filter With Position Aiding')
ylabel('error in degrees')
xlabel('time in minutes')
legend('lat err unaided', 'long err unaided', 'lat err KF', 'lon err KF',0)

figure
subplot(111)
plot(time/60,vel_l_err_unaided(:,1:2),time/60,vel_l_err_KF(:,1:2))
ylabel('velocity error in m/s')
xlabel('time in minutes')
legend('east unaided','north unaided', 'east KF', 'north KF',0)

figure
subplot(111)
plot(time/60,roll_err_unaided*180/pi,time/60,pitch_err_unaided*180/pi,...
     time/60,roll_err_KF*180/pi,time/60,pitch_err_KF*180/pi)
title('Euler Angle Errors')
ylabel('error in degrees')
xlabel('time in minutes')
legend('roll unaided','pitch unaided', 'roll KF', 'pitch KF',0)

figure
subplot(111)
plot(time/60,yaw_err_unaided*180/pi,time/60,yaw_err_KF*180/pi)
ylabel('yaw error in degrees')
xlabel('time in minutes')
legend('yaw unaided', 'yaw KF',0)

figure
plot(time/60,horz_pos_err_unaided,time/60,horz_pos_err_KF)
title('Horizontal Position Error')
ylabel('error in meters')
xlabel('time in minutes')
legend('unaided','KF',0)

figure
plot(time/60,horz_pos_err_KF)
title('Kalman Filter Horizontal Position Error')
ylabel('error in meters')
xlabel('time in minutes')

figure
plot(time/60,x_rms_KF,time/60,y_rms_KF)
title('Position Error Covariance Analysis')
ylabel('rms error in meters')
xlabel('time in minutes')
legend('x','y',0)

figure
plot(time/60,x_v_rms_KF,time/60,y_v_rms_KF)
title('Velocity Error Covariance Analysis')
ylabel('rms error in meters/sec')
xlabel('time in minutes')
legend('x','y',0)
