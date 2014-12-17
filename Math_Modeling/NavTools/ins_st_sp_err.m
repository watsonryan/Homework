%   
%  ins_st_sp_err.m
%
%  INS State Space Error modeling;  PSI-angle formulation
%
%  This m-file implements an INS error model state equation.
%
%  The continuous-time state equation is given by:
%
%     Xdot = F*X + G*U
%
%  where:  X is the state vector
%          Xdot is the state vector derivative
%          F is the continuous-time system dynamics matrix
%          U is the input vector
%          G is the connection matrix between the input vector and state
%          vector
%
%  There are nine state variables:
%
%  State vector is given by:
%  x(1) = x (east) position error
%  x(2) = y (north) position error
%  x(3) = z (up) position error
%  x(4) = x velocity error
%  x(5) = y velocity error
%  x(6) = z velocity error
%  x(7) = x component of psi-angle vector
%  x(8) = y component of psi-angle vector
%  x(9) = z component of psi-angle vector
%
%  There are six input parameters:
%
%  Input vector given by:
%  u(1) = body x (north) accel error
%  u(2) = body y (east) accel error
%  u(3) = body z (down) accel error
%  u(4) = body x (north) gyro error
%  u(5) = body y (east) gyro error
%  u(6) = body z (down) gyro error
%
%  This m-file is set up with the assumption of a stationary user.
%  However, with the trajectory generation capability of the INS
%  Toolbox,it could easily be extended to model a dynamic user.

%	References:
%   Huddle, J., "Inertial Navigation System Error Model Considerations
%               in Kalman Filtering Applications," in Control and Dynamic
%               Systems, Vol. 20, Academic Press, 1983.
%
%   Pue, A., "Integration of GPS with Inertial Navigation Systems,"
%             short course notes, Navtech Seminars, Springfield, VA, 2003.
%
%   Siouris, G., Aerospace Avionics Systems - A Modern Synthesis, Academic
%                Press, San Diego, 1993.
%
%  August 2004; August 2005
%  Copyright (c) 2004-2005 by GPSoft LLC
%  All Rights Reserved.
%
clear all
close all

T = 10;    % simulation time interval in seconds

accel_bias = (100e-6)*gravity(0,0);    % 100 micro-g bias
gyro_bias = (0.001)*(pi/180)*(1/3600);  % 0.001 deg/hr gyro bias converted to radians/sec

tru_lat = 45*pi/180;                   % set initial position
tru_long = 45*pi/180;                  % and wander angle
tru_alpha = 0;
tru_height = 0;
tru_vel_x = 0;
tru_vel_y = 0;

[rm,rp] = radicurv(tru_lat);  R = sqrt(rm*rp);   % earth radii of curvature


% In this section we specify the initial values of the state variables.
% Initial velocity and/or initial tilt errors can be specified here.
% Note that if initial position errors are modeled here the
% PSI-angle formulation couples position errors with attitude errors
%
% You can set any or all of the nine initial state variable values and/or
% the six input values to anything you wish.
     x = zeros(9,1);
     %x(5) = 0.1;  % initial north velocity error;
     %x(2) = (10/3600)*(pi/180)*R;  % 10 arc-second initial latitude error;
     %x(7) = (10/3600)*pi/180;  % 10 arc-second initial (east) tilt;
     %x(8) = -(10/3600)*pi/180;  % 10 arc-second initial (north) tilt;

     u = zeros(6,1);  
     %u(1) = accel_bias;  % body-x (north) accel bias;
     %u(2) = accel_bias;  % body-y (east) accel bias;
     %u(4) = gyro_bias;  % body-x gyro bias (north);
     u(5) = gyro_bias;  % body-y gyro bias (east);
     %u(6) = gyro_bias;  % body-z gyro bias (down);

time = 0:T:4*3600;    % 4-hour simulation
npts = length(time);

roll=0*pi/180;              % set initial roll to be zero radians
pitch=0*pi/180;            % set initial pitch
yaw=0*pi/180;              % set initial yaw

DCMnb=eulr2dcm([roll pitch yaw]);         % initialize nav-to-body
%                                        % direction cosine matrix
DCMbn = DCMnb';              % convert to body-to-nav DCM

lat_ex = tru_lat;
height_ex = tru_height;
vx_ex = tru_vel_x;
vy_ex = tru_vel_y;

C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

g = gravity(tru_lat,tru_height);

accel_vect_b = [0 0 -g]';     % true specific-force (acceleration) vector
accel_vect_L = C*DCMbn*accel_vect_b;


%                              % initialize earth-to-local-level
%                              % direction cosine matrix
DCMel = llw2dcm([tru_lat tru_long tru_alpha]);   

earthflg=0; vertmech=0;

%                              % calculate craft-rate in the L-frame
omega_el_L = crafrate(lat_ex,vx_ex,vy_ex,height_ex,DCMel,earthflg,vertmech);
%                              % convert craft-rate to the n-frame
omega_en_n = C*omega_el_L;

omega_ie_E = [0 0 7.292115e-5]';   % earth-rate expressed in the earth-frame
omega_ie_L = DCMel*omega_ie_E;   % earth-rate expressed in the L-frame
omega_ie_n = C*omega_ie_L;    % earth-rate expressed in the n-frame

%  Here we are forming the continuous-time system dynamics matrix
F11 = -1*skewsymm(omega_el_L);    % L-frame
F12 = eye(3);
F13 = 0*eye(3);
F21=eye(3); F21(1,1)=-g/R; F21(2,2)=-g/R; F21(3,3)=2*g/(R+tru_height);
F22 = -1*skewsymm(2*omega_ie_L + omega_el_L);   % L-frame
F23 = skewsymm(accel_vect_L);
F31 = 0*eye(3);
F32 = 0*eye(3);
F33 = (-1)*skewsymm(omega_ie_L + omega_el_L);   % L-frame
F = [F11 F12 F13; F21 F22 F23; F31 F32 F33];

%  Now form the input connection matrix, G
G11 = 0*eye(3);
G12 = 0*eye(3);
G21 = C*DCMbn;    % body-to-L_frame
G22 = 0*eye(3);
G31 = 0*eye(3);
G32 = (-1)*C*DCMbn;  % body-to-L_frame
G = [G11 G12; G21 G22; G31 G32];

[Fd,Gd] = cont2disc(F,G,T);

h = waitbar(0,'Looping over Time');
for i = 1:npts;
    x = Fd*x + Gd*u;
    x(6) = 0;  % emulate 'perfect' vertical aiding
    dr(i,:) = x(1:3)';
    dv(i,:) = x(4:6)';
    psi(i,:) = x(7:9)';
    
    waitbar(i/npts,h);
end
close(h)

theta(:,1) = -dr(:,2)/R;   % converting position error from linear
theta(:,2) = dr(:,1)/R;    % units to angular units
theta(:,3) = tan(tru_lat)*theta(:,2);  % this is for alpha=0 and ENU frame;

phi_angle = psi + theta;         % computing total attitude error

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % In the computation of 'est_DCMbn' below, the two multiplications
      % by 'C' deal with the fact that the phi_angle-vector has been computed
      % in the L-frame but our DCM is from body-to-nav
for i = 1:max(size(phi_angle)),
   est_DCMbn = C*(eye(3) - skewsymm(phi_angle(i,1:3)))*C*DCMbn;
   eul_vect = dcm2eulr(est_DCMbn);
   roll_err(i) = eul_vect(1);  % Note: Since the true values are equal
   pitch_err(i) = eul_vect(2); %  to zero, then the estimated values are
   yaw_err(i) = eul_vect(3);   %  equal to the errors
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:max(size(theta)),
    est_DCMel = (eye(3) - skewsymm(theta(i,1:3)))*DCMel;
    llw_vect = dcm2llw(est_DCMel);
    lat_err(i) = llw_vect(1) - tru_lat;
    lon_err(i) = llw_vect(2) - tru_long;
    alpha_err(i) = llw_vect(3) - tru_alpha;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_hours = time*(1/3600);

rad2arcsec = (180/pi)*3600;

figure
plot(t_hours,lat_err*rad2arcsec,t_hours,lon_err*rad2arcsec);
title('Inertial State Space Error Modeling')
ylabel('position error in arc-seconds')
xlabel('time in hours')
legend('lat error','lon error',0)

figure
plot(t_hours,theta(:,1)*rad2arcsec,t_hours,theta(:,2)*rad2arcsec);
title('Inertial State Space Error Modeling')
ylabel('theta error in arc-seconds')
xlabel('time in hours')
legend('thetax','thetay',0)

figure
plot(t_hours,dv(:,1),t_hours,dv(:,2))
title('Inertial State Space Error Modeling')
ylabel('velocity error in m/s')
xlabel('time in hours')
legend('dvx','dvy',0)

figure
subplot(211)
plot(t_hours,phi_angle(:,1)*rad2arcsec,t_hours,phi_angle(:,2)*rad2arcsec)
title('Inertial State Space Error Modeling')
ylabel('attitude error in arc-seconds')
legend('phix','phiy',0)
subplot(212)
plot(t_hours,phi_angle(:,3)*rad2arcsec)
legend('phiz',0)
xlabel('time in hours')

figure
subplot(211)
plot(t_hours,psi(:,1)*rad2arcsec,t_hours,psi(:,2)*rad2arcsec)
title('Inertial State Space Error Modeling')
ylabel('attitude error in arc-seconds')
legend('psix','psiy',0)
subplot(212)
plot(t_hours,psi(:,3)*rad2arcsec)
legend('psiz',0)
xlabel('time in hours')

figure
plot(t_hours,roll_err*rad2arcsec,t_hours,pitch_err*rad2arcsec,t_hours,yaw_err*rad2arcsec)
title('Inertial State Space Error Modeling')
ylabel('attitude error in arc-seconds')
xlabel('time in hours')
legend('roll','pitch','yaw',0)
