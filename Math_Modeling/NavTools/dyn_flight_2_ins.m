%
%   dyn_flight_2_ins.m
%
%   Same as DYN_FLIGHT_INS.M but here we are using a different
%   flight trajectory (see LOAD commands below)
%
%   Inertial Navigation System simulation with initialization and
%   sensor errors and a dynamic flight path
%
%   This program takes input from the flight trajectory generator
%   and a baseline INS run (see the 2 LOAD commands below).  The
%   baseline INS run (dyn_flt_ins_0.mat) is used to remove the
%   effects of imperfect INS algorithm integration (primarily due
%   to relatively low data rates) from the final error analysis.
%
%  September 2005
%  Copyright (c) 2005 by GPSoft LLC
%  All Rights Reserved.

clear all
close all

load dyn_flt_2_dat   % output from DYN_FLIGHT_2.M

load dyn_flt_2_ins_0   % output from DYN_FLIGHT_2_INS0.M

dph2rps = (pi/180)/3600;   % conversion constant from deg/hr to rad/sec

init_vel_e_err = 0.02;
init_vel_n_err = 0.02;   % 2 m/s initial north velocity error

%%init_x_tilt = (10/3600)*pi/180;  % 10 arc-second initial body-x tilt error
init_x_tilt = 0.0001;  % 0.1 milli-radian body-x tilt error
init_y_tilt = 0.0001;

g = gravity(0,0);
vxbias = 50.0e-6*g;       % 50 micro-g x-accel bias
vybias = 40.0e-6*g;       % 40 micro-g y-accel bias
vzbias = 0;

vxsferr = 0;        % 0 accel scale-factor errors
vysferr = 0;
vzsferr = 0;

vxstdev = 0;       % 0 accel white noise standard deviation
vystdev = 0;
vzstdev = 0;

thxbias = 0.01;       % 0 deg/hr x-gyro bias
thybias = 0.015;      % 0 deg/hr y-gyro bias
thzbias = 0.008;

thxsferr = 0;        % 0 gyro scale-factor errors
thysferr = 0;
thzsferr = 0;

thxstdev = 0;     % 0 deg/root-hour gyro white noise standard deviation
thystdev = 0;    % 0 deg/root-hour gyro white noise standard deviation
thzstdev = 0;

dvparam = [vxbias vybias vzbias;       % Set up parameter matrix for
   vxsferr vysferr vzsferr;            % GENDVERR
   vxstdev vystdev vzstdev];

dthparam = [thxbias thybias thzbias;   % Set up parameter matrix for
   thxsferr thysferr thzsferr;         % GENTHERR
   thxstdev thystdev thzstdev];

earthflg = 1;

orgllh = [39*pi/180 -82*pi/180 0];  % Start-point of ENU profile
[lat_prof,lon_prof,height_prof] = profconv(profile,orgllh);

npts = max(size(lat_prof));
vel_prof_L = profile(:,4:6);             % extracting velocity profile

DCMnb_prof = profile(:,10:18);           % extracting true nav-to-body
%                                        % direction cosine matrix

fprintf(1,' Generating Delta Thetas \n')
dthetbody = gendthet(DCMnb_prof);       % component of delta-theta associated
%                                      % with body motion

for k = 1:size(profile,1),
    dcmnb=[profile(k,10:12); profile(k,13:15); profile(k,16:18)];
    dcmbn=dcmnb';
    eulv=dcm2eulr(dcmbn);
    yaw_prof(k) = eulv(3);
end

fprintf(1,' Generating earth-to-local_level DCMs \n')
DCMel_prof = dcmelgen(lat_prof, lon_prof, yaw_prof);

deltaer = earthrot(time,DCMel_prof,DCMnb_prof);   % component of delta-theta
%                                                 % associated with earth-rate

deltacr = gendelcr2(lat_prof,vel_prof_L,height_prof,...
                       time,DCMnb_prof,DCMel_prof,earthflg);

deltheta = dthetbody + deltaer + deltacr;   % ideal (error-free) gyro output

dtherr = gentherr(deltheta,time,dthparam,98765);  % generate delta-theta errors

est_dtheta = deltheta + dtherr;   % form profiles of 'measured' delta-theta's

est_roll(1) = 0 + init_x_tilt;   % Aircraft is nominally level for the entire
est_pitch(1) = 0 + init_y_tilt;  % flight path.  Note that craft rate was
%                                % handled earlier.
est_yaw(1) = yaw_prof(1);

laterr=0;  longerr=0;  alphaerr=0;        % INITIALIZATION in this whole section
height = height_prof(1); height_err = 0;
height1 = height_prof(1);  height2 = height_prof(1);
est_lat(1) = lat_prof(1); est_lat(2) = lat_prof(1);
est_lon(1) = lon_prof(1);
vx1 = vel_prof_L(1,1) + init_vel_e_err; vx2 = vx1;
vy1 = vel_prof_L(1,2) + init_vel_n_err; vy2 = vy1;
vel_l(1,:) = [vx2 vy2 0];
vel2 = [vx1 vy1 0];  vel1 = vel2;
lat2 = lat_prof(1);  lat1 = lat_prof(1) - (lat_prof(2)-lat_prof(1));
est_DCMbn = ( eulr2dcm([est_roll(1) est_pitch(1) est_yaw(1)]) )';
est_DCMel = [DCMel_prof(1,1:3); DCMel_prof(1,4:6); DCMel_prof(1,7:9)];
vertmech = 0;
omega2_el_L = crafrate(lat_prof(1),vx1,vy1,height_prof(1),est_DCMel,earthflg,vertmech);

deltav_b = gendv(vel_prof_L,DCMnb_prof);    % Generate the component of delta-V
%                                           % associated with body motion relative
%                                           % to the earth

%                                           % Generate the Coriolis component
%                                           % of delta-V
dvcor = gendvcor2(lat_prof,vel_prof_L,height_prof,time,...
                          DCMnb_prof,DCMel_prof,earthflg);
dvtot = deltav_b + dvcor;

dverr = gendverr(dvtot,time,dvparam,76543);  % Generate delta-V errors

est_dv = dvtot + dverr;           % form profile of 'measured' delta-V's

fprintf(1,' . . . . . . . . . \n')
fprintf(1,' Starting nav computations \n')
C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU
est_DCMbn(:,:,1) = est_DCMbn;
est_DCMel(:,:,1) = est_DCMel;

h = waitbar(0,' Time Loop ');
for i = 2:npts,

   td12 = time(i) - time(i-1);
   tdex = 0.5*td12;
   tdint(i) = time(i) - time(i-1);
   
   est_DCMbn(:,:,i) = bodupdat(est_DCMbn(:,:,i-1),est_dtheta(i-1,1:3));
   
   [DCM_ll_I, omega_el_L(:,i), omega_ie_L(:,i)] = lclevupd(lat1,lat2,vx1,vx2,vy1,vy2,...
      height1,height2,td12,tdex,tdint(i),est_DCMel(:,:,i-1),vertmech,1,earthflg);
   
   est_DCMbn(:,:,i) = C*(DCM_ll_I*(C*est_DCMbn(:,:,i)));
   
   eul_vect = dcm2eulr(est_DCMbn(:,:,i));
   est_roll(i) = eul_vect(1);
   est_pitch(i) = eul_vect(2);
   est_yaw(i) = eul_vect(3);
   
   est_delv_b = est_dv(i-1,1:3);       % extract delta-V for current point in time
          
   del_Vl(:,i) = C*(est_DCMbn(:,:,i)*est_delv_b');

   omega1_el_L = omega2_el_L;   omega2_el_L = omega_el_L(:,i);
   [est_DCMel(:,:,i), DCM_ll_E] = navupdat(omega1_el_L,omega2_el_L,td12,est_DCMel(:,:,i-1),1);
   
   h_extr = extrapol(height1,height2,td12,tdex);
   lat_extr = extrapol(lat1,lat2,td12,tdex);
   g_extr(i) = gravity(lat_extr,h_extr);
 
   vtmp = velupdat(vel2,vel1,td12,tdex,del_Vl(:,i),...
      omega_el_L(:,i),est_DCMel(:,:,i),g_extr(i),0,tdint(i));
   
   vel_l(i,:) = vtmp';
   
   est_height(i,1) = height_prof(i);

   height1 = height2;  height2 = est_height(i,1);
    vx1 = vx2; vy1 = vy2;
    vx2 = vel_l(i,1);  vy2 = vel_l(i,2);
    vel1 = vel2;   vel2 = vel_l(i,:);
    llw_vect = dcm2llw(est_DCMel(:,:,i));
    est_lat(i) = llw_vect(1); est_lon(i) = llw_vect(2);  est_alpha(i) = llw_vect(3);
    lat1 = lat2;  lat2 = est_lat(i);
    est_lat_err(i) = est_lat(i)-est_lat0(i);
    est_lon_err(i) = est_lon(i)-est_lon0(i);
    
    waitbar(i/npts,h)
end
close(h)

% save dyn_flt_2_ins_dat lat_prof lon_prof height_prof vel_l est_roll est_pitch est_yaw ...
%     est_lat est_lon est_alpha del_Vl tdint omega_el_L omega_ie_L g_extr est_height ...
%     est_DCMbn est_DCMel est_lat_err est_lon_err npts time ...
%     est_alpha0 est_lat0 est_lon0 vel_l_0 est_roll0 est_pitch0 est_yaw0 -V6    % version 7.x

save dyn_flt_2_ins_dat lat_prof lon_prof height_prof vel_l est_roll est_pitch est_yaw ...
    est_lat est_lon est_alpha del_Vl tdint omega_el_L omega_ie_L g_extr est_height ...
    est_DCMbn est_DCMel est_lat_err est_lon_err npts time ...
    est_alpha0 est_lat0 est_lon0 vel_l_0 est_roll0 est_pitch0 est_yaw0    % version 6.x

