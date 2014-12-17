%
%   gps_ins_kf_fb.m
%
%   Inertial Navigation System simulation with initialization and
%   sensor errors and a dynamic flight path
%
%   After initialization, feedback the KF estimates of the accel
%   and gyro errors as corrections to subsequent measurements that
%   are then used in the INS pos/vel/att update
%
%  September 2009
%  Copyright (c) 2009 by GPSoft LLC
%  All Rights Reserved.

clear all
close all

load prof_1_ins_truth   % output from GEN_INS_TRUTH.M

dph2rps = (pi/180)/3600;   % conversion constant from deg/hr to rad/sec

earthflg = 1;      %  earth shape flag:  1 for ellipsoidal

C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

%%%%%%%% SPECIFICATION OF INS ERRORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

init_vel_e_err = 0.02;
init_vel_n_err = 0.02;   % 0.02 m/s initial north velocity error

init_x_tilt = 0.0001;  % 0.1 milli-radian body-x tilt error
init_y_tilt = 0.0001;
init_z_mis = 0*0.001;

g = gravity(0,0);
vxbias = 200e-6*g;       % 200 micro-g x-accel bias
vybias = 150e-6*g;       % 150 micro-g y-accel bias
vzbias = -100e-6*g;

vxsferr = 0;        % 0 accel scale-factor errors
vysferr = 0;
vzsferr = 0;

vxstdev = 0;       % 0 accel white noise standard deviation
vystdev = 0;
vzstdev = 0;

thxbias = 1;       % 1 deg/hr x-gyro bias
thybias = 1.5;      % 1.5 deg/hr y-gyro bias
thzbias = -2;

thxsferr = 0;        % 0 gyro scale-factor errors
thysferr = 0;
thzsferr = 0;

thxstdev = 0.09;     % 0.09 deg/root-hour gyro white noise standard deviation
thystdev = 0.05;    % 0.05 deg/root-hour gyro white noise standard deviation
thzstdev = -0.04;

dvparam = [vxbias vybias vzbias;       % Set up parameter matrix for
   vxsferr vysferr vzsferr;            % GENDVERR
   vxstdev vystdev vzstdev];

dthparam = [thxbias thybias thzbias;   % Set up parameter matrix for
   thxsferr thysferr thzsferr;         % GENTHERR
   thxstdev thystdev thzstdev];

%%%%%%%%%%%% END OF INS ERROR SPECIFICATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% Preallocate arrays for speed:  %%%%%%%%%%%%%%%%%%%%%%%%%%
accum_gyro_x = zeros(npts,1);
accum_gyro_y = zeros(npts,1);
accum_gyro_z = zeros(npts,1);
accum_accel_x = zeros(npts,1);
accum_accel_y = zeros(npts,1);
accum_accel_z = zeros(npts,1);
est_lat = zeros(npts,1);
est_lon = zeros(npts,1);
est_alpha = zeros(npts,1);
est_height = zeros(npts,1);
est_roll = zeros(npts,1);
est_pitch = zeros(npts,1);
est_yaw = zeros(npts,1);

state_est = zeros(18,npts);
x_rms_KF = zeros(npts,1);
y_rms_KF = zeros(npts,1);
z_rms_KF = zeros(npts,1);
x_v_rms_KF = zeros(npts,1);
y_v_rms_KF = zeros(npts,1);
z_v_rms_KF = zeros(npts,1);
x_psi_rms_KF = zeros(npts,1);
y_psi_rms_KF = zeros(npts,1);
z_psi_rms_KF = zeros(npts,1);
x_accel_rms_KF = zeros(npts,1);
y_accel_rms_KF = zeros(npts,1);
z_accel_rms_KF = zeros(npts,1);
x_gyro_rms_KF = zeros(npts,1);
y_gyro_rms_KF = zeros(npts,1);
z_gyro_rms_KF = zeros(npts,1);
est_roll_KF = zeros(npts,1);
est_pitch_KF = zeros(npts,1);
est_yaw_KF = zeros(npts,1);
est_lat_KF = zeros(npts,1);
est_lon_KF = zeros(npts,1);
est_alpha_KF = zeros(npts,1);
est_height_KF = zeros(npts,1);
est_vel_l_KF =  zeros(npts,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% INS MECHANIZATION ALGORITHM INITIALIZATION %%%%%%%%%%%%%%%%%

est_roll(1) = tru_roll(1) + init_x_tilt;
est_pitch(1) = tru_pitch(2) + init_y_tilt;
est_yaw(1) = tru_yaw(3) + init_z_mis;
    
laterr=0;  longerr=0;  alphaerr=0; 
est_alpha(1) = tru_alpha(1);
est_height(1) = tru_height(1);
height1 = tru_height(1);  height2 = tru_height(1);
est_lat(1) = tru_lat(1); est_lat(2) = tru_lat(1);
est_lon(1) = tru_lon(1);
vx1 = tru_vel_L(1,1) + init_vel_e_err; vx2 = vx1;
vy1 = tru_vel_L(1,2) + init_vel_n_err; vy2 = vy1;
vel_l(1,:) = [vx2 vy2 0];
vel2 = [vx1 vy1 0];  vel1 = vel2;
lat2 = est_lat(1);  lat1 = est_lat(1) - (est_lat(2)-est_lat(1));
est_DCMbn = ( eulr2dcm([est_roll(1) est_pitch(1) est_yaw(1)]) )';
est_DCMel = llw2dcm([est_lat(1) est_lon(1) est_alpha(1)]);
vertmech = 0;
omega2_el_L = crafrate(est_lat(1),vx1,vy1,est_height(1),est_DCMel,earthflg,vertmech);

%%%%%%%%%%%% END INS MECHANIZATION INITIALIZATION  %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% BEGIN GENERATION OF IMU MEASUREMENTS %%%%%%%%%%%%%%%%%%%%%%%

dtherr = gentherr(deltheta,time,dthparam,98765);  % generate delta-theta errors

est_dtheta = deltheta + dtherr;   % form profiles of 'measured' delta-theta's

dverr = gendverr(dvtot,time,dvparam,76543);  % Generate delta-V errors

est_dv = dvtot + dverr;           % form profile of 'measured' delta-V's

%%%%%%%%%%% END GENERATION OF IMU MEASUREMENTS %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% BEGIN GPS INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mpmat=mpgen(24,3600,1,54321);  % simulate GPS multipath error
loadgps  % load the GPS constellation parameters into global memory
%%gpstime = 40000;  % specify the GPS time-of-day for the simulation
gpstime = 22000;  % specify the GPS time-of-day for the simulation
%%%%%%%%%%% END GPS INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% INSERT KF INITIALIZATION PARAMETERS %%%%%%%%%%%%%%%%

H = zeros(3,18);  H(1,1) = 1;  H(2,2) = 1;  H(3,3) = 1;
%                                 % The measurement vector z consists of
%                                 % an x, y and z position differences only.
%                                 % Since we are working in the L frame
%                                 % this corresponds to east, north and up
%                                 % position differences (since alpha=0).

%  Defining the components of the initial prediction error covariance
%  matrix
%%init_pos_rms = 3;  % initial position error (RMS in meters)
init_pos_rms = 1.5;  % initial position error (RMS in meters)
Peast_pos = init_pos_rms^2;  Pnorth_pos = init_pos_rms^2;
                Pup_pos = init_pos_rms^2;
Peast_vel = 2^2;  Pnorth_vel = 2^2;  Pup_vel = 2^2;
Ppsi_x = 0.0001^2;  Ppsi_y = 0.0001^2;  Ppsi_z = 0.0001^2;
Pacc_x = (100*9.81e-6)^2; 
Pacc_y = (100*9.81e-6)^2; 
Pacc_z = (100*9.81e-6)^2;
Pgyr_x = (0.05*dph2rps)^2; 
Pgyr_y = (0.05*dph2rps)^2; 
Pgyr_z = (0.05*dph2rps)^2;

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
update = 0;  count = 0;

X_pre = zeros(18,1);   % Initialize the state prediction

est_lat_KF(1) = est_lat(1);
est_lon_KF(1) = est_lon(1);
est_height_KF(1) = est_height(1);
est_alpha_KF(1) = est_alpha(1);
est_vel_l_KF(1,:) =  vel_l(1,:);
est_roll_KF(1) = est_roll(1);
est_pitch_KF(1) = est_pitch(1);
est_yaw_KF(1) = est_yaw(1);

vertmech = 0;
earthflg = 1;
omega_ie_E = [0 0 7.292115e-5]';
init_corrections = 0;
prev_KF_update = 0;
KF_corrections = 0;

%%%%%%%%%%%%%%%%%%%% END OF KALMAN FILTER INITIALIZATION %%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%% BEGIN TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

accum_accel_bias = zeros(3,1);
accum_gyro_bias = zeros(3,1);

fprintf(1,' . . . . . . . . . \n')
fprintf(1,' Starting nav computations \n')

npts = length(time);
h = waitbar(0,' Time Loop ');
for i = 2:npts,

   td12 = time(i) - time(i-1);
   tdex = 0.5*td12;
   tdint = time(i) - time(i-1);
   
   est_dtheta(i-1,1:3) = est_dtheta(i-1,1:3) - (accum_gyro_bias*td12)';
        accum_gyro_x(i) = accum_gyro_bias(1);
        accum_gyro_y(i) = accum_gyro_bias(2);
        accum_gyro_z(i) = accum_gyro_bias(3);
   
   est_DCMbn = bodupdat(est_DCMbn,est_dtheta(i-1,1:3));
   
   [DCM_ll_I, omega_el_L, omega_ie_L] = lclevupd(lat1,lat2,vx1,vx2,vy1,vy2,...
      height1,height2,td12,tdex,tdint,est_DCMel,vertmech,1,earthflg);
   
   est_DCMbn = C*(DCM_ll_I*(C*est_DCMbn));
   
   eul_vect = dcm2eulr(est_DCMbn);
   est_roll(i) = eul_vect(1);
   est_pitch(i) = eul_vect(2);
   est_yaw(i) = eul_vect(3);
   
   est_delv_b = est_dv(i-1,1:3) - (accum_accel_bias*td12)';
        accum_accel_x(i) = accum_accel_bias(1);
        accum_accel_y(i) = accum_accel_bias(2);
        accum_accel_z(i) = accum_accel_bias(3);
          
   del_Vl = C*(est_DCMbn*est_delv_b');

   omega1_el_L = omega2_el_L;   omega2_el_L = omega_el_L;
   [est_DCMel, DCM_ll_E] = navupdat(omega1_el_L,omega2_el_L,td12,est_DCMel,1);
   
   h_extr = extrapol(height1,height2,td12,tdex);
   lat_extr = extrapol(lat1,lat2,td12,tdex);
   g_extr = gravity(lat_extr,h_extr);
 
   vtmp = velupdat(vel2,vel1,td12,tdex,del_Vl,...
      omega_el_L,est_DCMel,g_extr,0,tdint);
   
   vel_l(i,:) = vtmp';
   
   est_height(i,1)= est_height(i-1,1) + tdint*mean([vel_l(i,3); vel_l(i-1,3)]);

   height1 = height2;  height2 = est_height(i,1);
    vx1 = vx2; vy1 = vy2;
    vx2 = vel_l(i,1);  vy2 = vel_l(i,2);
    vel1 = vel2;   vel2 = vel_l(i,:);
    llw_vect = dcm2llw(est_DCMel);
    est_lat(i) = llw_vect(1); est_lon(i) = llw_vect(2);  est_alpha(i) = llw_vect(3);
    lat1 = lat2;  lat2 = est_lat(i);
    
    %%%%%%%%%%%%%% BEGIN KALMAN FILTER ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%
    
    % Forming the system dynamics matrix
    [rm,rp] = radicurv(est_lat(i));  radius_e = sqrt(rm*rp);
    accel_vect_L = del_Vl*(1/tdint);
    est_DCMel = llw2dcm([est_lat(i) est_lon(i) est_alpha(i)]);
    omega_el_L = crafrate(est_lat(i),vel_l(i,1),vel_l(i,2),...
        est_height(i),est_DCMel,earthflg,vertmech);
    F11 = -1*skewsymm(omega_el_L);    % L-frame
    F12 = eye(3);
    F13 = 0*eye(3);
    grav = gravity(est_lat(i),est_height(i));
    F21=eye(3); F21(1,1)=-grav/radius_e; 
    F21(2,2)=-grav/radius_e; F21(3,3)=2*grav/(radius_e+est_height(i,1));
    omega_ie_L = est_DCMel*omega_ie_E;
    F22 = -1*skewsymm(2*omega_ie_L + omega_el_L);   % L-frame
    F23 = skewsymm(accel_vect_L);
    F31 = 0*eye(3);
    F32 = 0*eye(3);
    F33 = (-1)*skewsymm(omega_ie_L + omega_el_L);   % L-frame
    F_ins = [F11 F12 F13; F21 F22 F23; F31 F32 F33];

    F = zeros(18,18); F(1:9,1:9) = F_ins; 
    est_DCMbn = ( eulr2dcm([est_roll(i) est_pitch(i) est_yaw(i)]) )';
    F(10:12,10:12) = (-1/tau_accel)*eye(3);  F(4:6,10:12) = C*est_DCMbn;
    F(13:15,13:15) = (-1/tau_gyro)*eye(3);  F(7:9,13:15) = -1*C*est_DCMbn;

    
    %  Numerical calculation of the state transition matrix (PHI) and
    %  system noise covariance matrix (Q)
    A = zeros(36,36);
    A(1:18,1:18) = -1*F;
    A(1:18,19:36) = G*W*G';
    A(19:36,19:36) = F';
    A = A*tdint;

    B = expm(A);

    PHI_trans = B(19:36,19:36);
    PHI = PHI_trans';

    Q = PHI*B(1:18,19:36);

    % Calculate the Kalman gain
    % THE NEXT LINE IMPLEMENTS THE TRADITIONAL KALMAN GAIN:
    K = (P_pre_KF*H')/(H*P_pre_KF*H' + R);
    
    % keep track of elapsed time (in seconds) since last Kalman update
   count = count + tdint;

    if count >= .1,    %Hz updates
        update = 1;
        update_interval = count;
        count = 0;
    end
    
    if update == 1,
        %%  Simulate GPS Receiver output here
        gpstime = gpstime + update_interval;
        tru_pos_xyz = llh2xyz([tru_lat(i) tru_lon(i) tru_height(i)]);
        [svxyzmat,svid] = gensv(tru_pos_xyz,gpstime,5);
        %%[prvec,adrvec] = genrng(1,tru_pos_xyz,svxyzmat,svid,gpstime,[1 1 0 1 1],[],mpmat);
        [prvec,adrvec] = genrng(1,tru_pos_xyz,svxyzmat,svid,gpstime,[1 0 0 0 0],[],mpmat);
        gps_xyzt = olspos(prvec,svxyzmat);
        gps_llh = xyz2llh(gps_xyzt(1:3));

        ins_lat = est_lat(i);
        ins_lon = est_lon(i);
        ins_height = est_height(i);

        %  Kalman observation vector (Z) formed by position differences.
        %  Angular position differences converted to linear differences.
        Z(1,1) = (ins_lon - gps_llh(2))*rp*cos(ins_lat);  % east component
        Z(2,1) = (ins_lat - gps_llh(1))*rm;               % north component
        Z(3,1) = ins_height - gps_llh(3);             % up component

        X_est = X_pre + K*(Z - H*X_pre);

        P_est_KF = (eye(18) - K*H)*P_pre_KF;

        prev_KF_update = 1;
        update = 0;
    else
        X_est = zeros(18,1);
        P_est_KF = P_pre_KF;
    end
    X_pre = PHI*X_est;      
    P_pre_KF = PHI*P_est_KF*PHI' + Q;
    
    %%%%%%%%%%% SAVE VARIABLES FOR PLOTTING LATER %%%%%%%%%%%%%%%%%%%%%
    state_est(:,i) = X_est;
    x_rms_KF(i) = sqrt( P_est_KF(1,1) );
    y_rms_KF(i) = sqrt( P_est_KF(2,2) );
    z_rms_KF(i) = sqrt( P_est_KF(3,3) );
    x_v_rms_KF(i) = sqrt( P_est_KF(4,4) );
    y_v_rms_KF(i) = sqrt( P_est_KF(5,5) );
    z_v_rms_KF(i) = sqrt( P_est_KF(6,6) );
    x_psi_rms_KF(i) = sqrt( P_est_KF(7,7) );
    y_psi_rms_KF(i) = sqrt( P_est_KF(8,8) );
    z_psi_rms_KF(i) = sqrt( P_est_KF(9,9) );
    x_accel_rms_KF(i) = sqrt( P_est_KF(10,10) );
    y_accel_rms_KF(i) = sqrt( P_est_KF(11,11) );
    z_accel_rms_KF(i) = sqrt( P_est_KF(12,12) );
    x_gyro_rms_KF(i) = sqrt( P_est_KF(13,13) );
    y_gyro_rms_KF(i) = sqrt( P_est_KF(14,14) );
    z_gyro_rms_KF(i) = sqrt( P_est_KF(15,15) );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta(1,1) = -X_est(2)/rm;   % converting position error from linear
    theta(2,1) = X_est(1)/rp;    % units to angular units
    theta(3,1) = tan(est_lat(i))*theta(2);  % this is for alpha=0 and ENU frame;

    psi = X_est(7:9);
    phi_angle = psi + theta;         % computing total attitude error

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % This is the Kalman estimate of the true DCM
    est_DCMbn_KF = C*(eye(3) + skewsymm(phi_angle))*C*est_DCMbn; 
    
    eul_vect = dcm2eulr(est_DCMbn_KF);
    est_roll_KF(i) = eul_vect(1);
    est_pitch_KF(i) = eul_vect(2);
    est_yaw_KF(i) = eul_vect(3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % This is the Kalman estimate of the true DCM
    est_DCMel_KF = (eye(3) + skewsymm(theta))*est_DCMel;
    
    llw_vect = dcm2llw(est_DCMel_KF);
    est_lat_KF(i) = llw_vect(1);
    est_lon_KF(i) = llw_vect(2);
    est_alpha_KF(i) = llw_vect(3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    est_height_KF(i,1) = est_height(i) - X_est(3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    est_vel_l_KF(i,:) = ( (eye(3) + skewsymm(theta))*(vel_l(i,:)' - X_est(4:6)) )';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%% END KALMAN FILTER ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
   %%%%%%% THE FOLLOWING SECTION INCORPORATES THE KALMAN FILTER CORRECTIONS
   %%%%%%% INTO THE INS POS/VEL/ATT UPDATE ALGORITHM 
   if (prev_KF_update==1), 
       time_in_minutes = time(i)/60;
       accum_accel_bias = accum_accel_bias + X_est(10:12);
       accum_gyro_bias = accum_gyro_bias + X_est(13:15);
       est_DCMbn = est_DCMbn_KF;
       vel_l(i,:) = est_vel_l_KF(i,:);
       vel_l(i-1,:) = est_vel_l_KF(i-1,:);
       est_DCMel = est_DCMel_KF;
       est_height(i) = est_height_KF(i,1);
       est_height(i-1) = est_height_KF(i-1,1);

       vel1 = vel_l(i,:);
       vel2 = vel_l(i,:);
       vx1 = vel_l(i,1); vy1 = vel_l(i,2);
       vx2 = vel_l(i,1); vy2 = vel_l(i,2);
       lat1 = est_lat_KF(i);
       lat2 = est_lat_KF(i);
       height2 = est_height(i);
       height1 = est_height(i);

       X_pre = zeros(18,1);
       KF_corrections = 1;
       prev_KF_update = 0;
   end
   %%%%%%%%%%%%% END OF KALMAN FILTER CORRECTION FEEDBACK %%%%%%%%%%%%%%

   if rem(i,1000)==0,
      waitbar(i/npts,h)
   end
end
close(h)

%%%%%%%%%%%%%%%% END OF TIME LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h = waitbar(0,' Computing Position, Velocity, Attitude Errors ');
east_pos_err_KF = zeros(npts,1);
north_pos_err_KF = zeros(npts,1);
up_pos_err_KF = zeros(npts,1);
for i = 1:npts,
   truxyz = llh2xyz([tru_lat(i) tru_lon(i) tru_height(i)]); 
   insxyz_KF = llh2xyz([est_lat_KF(i) est_lon_KF(i) est_height_KF(i)]);
   enu = xyz2enu(insxyz_KF,truxyz);
   east_pos_err_KF(i) = enu(1);
   north_pos_err_KF(i) = enu(2);
   up_pos_err_KF(i) = enu(3);

   if rem(i,1000)==0,
      waitbar(i/npts,h)
   end
end
close(h)

lat_err_KF = est_lat_KF-tru_lat;
lon_err_KF = est_lon_KF-tru_lon;
alpha_err_KF = est_alpha_KF - tru_alpha;
vel_l_err_KF = est_vel_l_KF - tru_vel_L;
roll_err_KF = est_roll_KF - tru_roll;
pitch_err_KF = est_pitch_KF - tru_pitch;
yaw_err_KF = est_yaw_KF - tru_yaw;

close all

time_end = max(time/60);

save kf_fb_out

plot_kf_fb

