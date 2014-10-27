%
%   insdem05.m
%      Static user with north accel bias
%
%      This is effectively a single-channel simulation
%      to observe Schuler oscillations
%
%      - Earth-rate set to zero in local-level update
%      - Coriolis correction NOT included in velocity update
%

clear all
close all

%n_accel_bias = (4.848e-5)*gravity(0,0);    % equivalent to 10 arc-sec tilt
n_accel_bias = (100e-6)*gravity(0,0);    % 100 micro-g bias
deltat = 10;                         % simulation time step: 10 seconds
time = 0:deltat:4*3600;
npts = max(size(time));

phi=0*pi/180;              % set initial roll to be zero radians
theta=0*pi/180;            % set initial pitch
psi=0*pi/180;              % set initial yaw

DCMnb=eulr2dcm([phi theta psi]);         % initialize nav-to-body
%                                        % direction cosine matrix
DCMbn = DCMnb';              % convert to body-to-nav DCM

est_roll(1) = phi;            % set initial Euler angle
est_pitch(1) = theta;         % estimates 
est_yaw(1) = psi;

tru_lat = 45*pi/180;                   % set initial position
tru_long = 45*pi/180;                  % and wander angle
tru_alpha = 0;

%                              % initialize earth-to-local-level
%                              % direction cosine matrix
DCMel = llw2dcm([tru_lat tru_long tru_alpha]);   

est_lat = tru_lat;             % set initial position and
est_long = tru_long;           % wander angle estimates
est_alpha = tru_alpha;

lat1 = est_lat;              % initialize current (2) and previous (1)
lat2 = est_lat;              % values of latitude (used in LCLEVUPD)

vel_L(1,:) = [0 0 0];         % initialize velocity vector estimate in 
%                             % the local-level frame

vel2 = vel_L(1,:);  vel1 = vel2;   % initialize current and previous
%                                  % velocity vectors (used in LCLEVUPD)

vx1 = vel_L(1,1);         % initialize current and previous values
vx2 = vel_L(1,1);          % of east (x) and north (y) velocity components
vy1 = vel_L(1,2);         % (this demo uses a north-slaved mechanization)
vy2 = vel_L(1,2);

height1 = 0; height2 = 0;    % initialize current and previous values
%                            % of altitude

vertmech = 0;          % specify vertical craft rate mechanization
%                      % (0 is north-slaved)

earthflg = 0;            % specify earth shape flag (0 = spherical;
%                        % 1 = ellipsoidal)

n_bias = n_accel_bias*deltat;      % form delta-V component due
%                                  % to accelerometer bias

g = gravity(tru_lat,height1); % determine value of plumb-bob gravity
%                             % at the vehicle location

accum_g = g*deltat;           % form delta-V component due to gravity

delv_b = [n_bias 0 -accum_g];   % form delta-V vector in body frame

C = [0 1 0; 1 0 0; 0 0 -1];    % conversion matrix between NED and ENU

for i = 2:npts,
   td12 = time(i) - time(i-1);
   tdex = 0.5*td12;
   tdint = td12;
   %                           % Compute the direction cosine matrix
   %                           % relating the angular displacement of
   %                           % of the local-level frame over the
   %                           % update interval
   %                               % Special version of LCLEVUPD used here
   %                               % since we are not yet simulating 
   %                               % the effects of earth rate
 	[DCM_ll_I,omega_el_L] = lclevtmp(lat1,lat2,vx1,vx2,vy1,vy2,...
                  height1,height2,td12,tdex,tdint,DCMel,vertmech,1,earthflg);

   DCMbn = C*(DCM_ll_I*(C*DCMbn));        % Update the body-to-nav DCM to
   %                                      % account for the motion of the
   %                                      % local-level frame
   
   eul_vect = dcm2eulr(DCMbn);      % compute Euler angle vector
   est_roll(i) = eul_vect(1);
   est_pitch(i) = eul_vect(2);
   est_yaw(i) = eul_vect(3);
   
   del_VL = C*(DCMbn*delv_b');     % Convert delta-V vector from the body
   %                               % frame to the local-level frame
   
   if i == 2,
     omega1_el_L = omega_el_L;
   else
      omega1_el_L = omega2_el_L;     % Save current and previous craft
   end                               % rate vectors for NAVUPDAT
   omega2_el_L = omega_el_L;
   
   %                                % Update earth-to-local-level DCM
   [DCMel, DCM_ll_E] = navupdat(omega1_el_L,omega2_el_L,td12,DCMel,1);
   %%[DCMel, DCM_ll_E] = navupd2(omega1_el_L,omega2_el_L,td12,DCMel,1);

   %                             % Update velocity
   vtmp = velupdat(vel2,vel1,td12,tdex,del_VL,omega_el_L,DCMel,g,1,td12);

   vel_L(i,:) = vtmp';         % save velocity estimate
   
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
plot(time/3600,(est_lat-tru_lat)*3600*180/pi)
title('INSDEM05 - Static User at 45 deg Lat - 100 micro-g North Accel Bias')
ylabel('latitude error in arc-seconds')
xlabel('time in hours')

subplot(212)
plot(time/3600,vel_L(:,2))
ylabel('north velocity error in m/s')
xlabel('time in hours')

figure
plot(time/3600,est_roll*3600*180/pi,...
    time/3600,est_pitch*3600*180/pi,...
    time/3600,est_yaw*3600*180/pi)
ylabel('attitude errors in arc-seconds')
xlabel('time in hours')
legend('roll','pitch','yaw',0)
