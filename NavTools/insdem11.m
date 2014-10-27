%
%   insdem11.m
%      Static user with east gyro noise
%
%      - All other errors are set to zero
%

clear all
close all

dprh2rprs = (pi/180)/sqrt(3600);  % conversion factor going from
%                                 % degrees-per-root-hour to
%                                 % radians-per-root-second

%                  % 0.01 deg/square-root-of-hour white noise standard deviation
n_gyro_noise_stdev = 0*dprh2rprs;
e_gyro_noise_stdev = 0.01*dprh2rprs;
d_gyro_noise_stdev = 0*dprh2rprs;

n_tilt = 0;
e_tilt = 0;

n_init_vel_error = 0;
e_init_vel_error = 0;

dph2rps = (pi/180)/3600;  % conversion factor going from
%                         % degrees-per-hour to radians-per-second

n_gyro_bias = 0*dph2rps;
e_gyro_bias = 0*dph2rps;
d_gyro_bias = 0*dph2rps;

n_accel_bias = 0;
e_accel_bias = 0;
d_accel_bias = 0;

deltat = 10;
time = 0:deltat:4*3600;
npts = max(size(time));

randn('seed',54321)              % Set random number generator seed
%                                 % for repeatibility
%
%                                 % Generate delta-theta component
%                                 % associated with gyro noise
n_gyro_noise = n_gyro_noise_stdev*sqrt(deltat)*randn(npts,1);
e_gyro_noise = e_gyro_noise_stdev*sqrt(deltat)*randn(npts,1);
d_gyro_noise = d_gyro_noise_stdev*sqrt(deltat)*randn(npts,1);
gyro_noise = [n_gyro_noise e_gyro_noise d_gyro_noise];

phi=0*pi/180 + n_tilt;
theta=0*pi/180 + e_tilt;
psi=0*pi/180;
DCMnb=eulr2dcm([phi theta psi]);  DCMbn = DCMnb';

roll(1) = phi; pitch(1) = theta; yaw(1) = psi;

tru_lat = 45*pi/180;
tru_long = 45*pi/180;
tru_alpha = 0;
DCMel = llw2dcm([tru_lat tru_long tru_alpha]);

est_lat = tru_lat;
est_long = tru_long;
est_alpha = tru_alpha;

lat1 = tru_lat; 
lat2 = tru_lat;

vx1 = 0+e_init_vel_error;             % Incorporate initial velocity
vx2 = vx1;                            % errors into initialization
vy1= 0+n_init_vel_error;              % parameters
vy2 = vy1;
vel_L(1,:) = [vx2 vy2 0];
vel2 = vel_L(1,:);  vel1 = vel2;

height1 = 0; height2 = 0;
vertmech = 0;
earthflg = 1;

n_bias = n_accel_bias*deltat;
e_bias = e_accel_bias*deltat;
d_bias = d_accel_bias*deltat;
g = gravity(0,0);
accum_g = g*deltat;
delv_b = [n_bias e_bias (-accum_g+d_bias)];

DCMel_prof(1,1:3) = DCMel(1,1:3);              % Very simple profile since the
DCMel_prof(1,4:6) = DCMel(2,1:3);              % true user attitude remains at
DCMel_prof(1,7:9) = DCMel(3,1:3);              % [0 0 0] and the true user
DCMel_prof = [DCMel_prof; DCMel_prof];         % position is stationary at
DCMnb_prof(1,1:3) = DCMnb(1,1:3);              % [0 0 0]
DCMnb_prof(1,4:6) = DCMnb(2,1:3);
DCMnb_prof(1,7:9) = DCMnb(3,1:3);
DCMnb_prof = [DCMnb_prof; DCMnb_prof];

deltaer = earthrot([deltat; 2*deltat],DCMel_prof,DCMnb_prof);

acc_n_gyro_bias = n_gyro_bias*deltat;
acc_e_gyro_bias = e_gyro_bias*deltat;
acc_d_gyro_bias = d_gyro_bias*deltat; 

dtherr = [acc_n_gyro_bias acc_e_gyro_bias acc_d_gyro_bias];

est_dtheta = deltaer + dtherr;     % Form the 'measured' delta-theta
%                                  % vector.  Since we are simulating
%                                  % a static vehicle, the delta-theta
%                                  % components are constant EXCEPT NOISE.
%                                  % Accordingly, the noise has to be
%                                  % added inside the time loop.

C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

for i = 2:npts,
   td12 = time(i) - time(i-1);
   tdex = 0.5*td12;
   tdint = td12;
   %                           % add current value of gyro noise to the
   %                           % delta-theta measurement
   new_est_dtheta = est_dtheta + gyro_noise(i,:);
   
   DCMbn = bodupdat(DCMbn,new_est_dtheta);   % Update the body-to-nav DCM with
   %                                         % the gyro output.
   
   %                           % Compute the direction cosine matrix
   %                           % relating the angular displacement of
   %                           % the local level frame over the update
   %                           % interval.  Note that this accounts for
   %                           % earth-rate and craft-rate.
 	[DCM_ll_I,omega_el_L,omega_ie_L] = lclevupd(lat1,lat2,vx1,vx2,vy1,vy2,...
             height1,height2,td12,tdex,tdint,DCMel,vertmech,1,earthflg);
   %
   DCMbn = C*(DCM_ll_I*(C*DCMbn));     % Update the body-to-nav DCM to
   %                                % account for the motion of the
   %                                % local-level frame
   
   eul_vect = dcm2eulr(DCMbn);
   roll(i) = eul_vect(1);
   pitch(i) = eul_vect(2);
   yaw(i) = eul_vect(3);
   
   del_VL = C*(DCMbn*delv_b');    % Convert delta-V vector from the body
   %                              % frame to the local-level frame

   if i == 2,
     omega1_el_L = omega_el_L;
   else
      omega1_el_L = omega2_el_L;
   end
   omega2_el_L = omega_el_L;
   [DCMel, DCM_ll_E] = navupdat(omega1_el_L,omega2_el_L,td12,DCMel,1);
   %%[DCMel, DCM_ll_E] = navupd2(omega1_el_L,omega2_el_L,td12,DCMel,1);

   %                              % Update velocity
   vtmp = velupdat(vel2,vel1,td12,tdex,del_VL,omega_el_L,DCMel,g,0,td12);

   vel_L(i,:) = vtmp';
   
   vx1 = vx2; vy1 = vy2;
    vx2 = vel_L(i,1);  vy2 = vel_L(i,2);
    vel1 = vel2;  vel2 = vel_L(i,:);
    llw_vect = dcm2llw(DCMel);  lat1 = lat2;  lat2 = llw_vect(1);
    est_lat(i) = llw_vect(1); 
    est_long(i) = llw_vect(2);  
    est_alpha(i) = llw_vect(3);
end

close all
subplot(211)
lat_err = (est_lat-tru_lat)*(180/pi)*3600;
long_err = (est_long-tru_long)*(180/pi)*3600;
plot(time/3600,lat_err,time/3600,long_err)
title('INSDEM11: 0.01 deg/root-hr East Gyro Noise')
ylabel('error in arc-seconds')
xlabel('time in hours')
text(1.2,-75,'LAT')
text(2.1,15,'LONG')
grid

subplot(212)
plot(time/3600,vel_L(:,1),time/3600,vel_L(:,2))
ylabel('velocity error in m/s')
xlabel('time in hours')
text(3.6,0.5,'EAST')
text(2.05,1.5,'NORTH')
grid

figure
plot(time/3600,roll*3600*180/pi,time/3600,pitch*3600*180/pi)
title('INSDEM11: 0.01 deg/root-hr East Gyro Noise')
ylabel('roll and pitch error in arc-seconds')
xlabel('time in hours')
text(1.8,-60,'pitch')
text(2.8,20,'roll')

figure
plot(time/3600,yaw*3600*180/pi)
title('INSDEM11: 0.01 deg/root-hr East Gyro Noise')
ylabel('yaw error in arc-seconds')
xlabel('time in hours')

