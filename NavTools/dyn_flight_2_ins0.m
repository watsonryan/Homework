%
%   dyn_flight_2_ins0.m
%
%   Same as DYN_FLIGHT_INS0.M but here we are using a different
%   flight trajectory (see LOAD command below)
%
%   Baseline (i.e., no INS sensor or initialization errors) run 
%
%   PURPOSE OF THIS PROGRAM:  Even if no initialization or sensor errors
%   are modeled, the INS updating algorithms still are not perfect.  This
%   is particularly true when relatively low data rates are employed in a
%   simulation (i.e., 1 Hz).  By running a baseline without initialization
%   and sensor errors, the resulting pos/vel/att error can be removed later
%   when performing a full simulation.  This way the impact of
%   initialization and sensor errors can be isolated from the underlying
%   algorithmic integration errors.
%
%   - INS position/velocity/attitude updating
%   - Simplified F-16 model dynamics
%   - Local-level trajectory is converted to earth-referenced
%     trajectory prior to generation/simulation of gyro and
%     accel outputs and subsequent INS processing
%
%
%  August 2005
%  Copyright (c) 2005 by GPSoft LLC
%  All Rights Reserved.
%

clear all
close all

load dyn_flt_2_dat    % output from DYN_FLIGHT_2.M

dph2rps = (pi/180)/3600;   % conversion constant from deg/hr to rad/sec

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

yaw_prof = yaw*pi/180;

fprintf(1,' Generating earth-to-local_level DCMs \n')
DCMel_prof = dcmelgen(lat_prof, lon_prof, yaw_prof);

deltaer = earthrot(time,DCMel_prof,DCMnb_prof);   % component of delta-theta
%                                                 % associated with earth-rate

%                                    % generate the component of delta-theta
%                                    % associated with craft rate
deltacr = gendelcr2(lat_prof,vel_prof_L,height_prof,...
                       time,DCMnb_prof,DCMel_prof,earthflg);

deltheta = dthetbody + deltaer + deltacr;   % ideal (error-free) gyro output

dtherr = zeros(size(deltheta));        % no delta-theta errors for this demo

est_dtheta = deltheta + dtherr;   % form profiles of 'measured' delta-theta's

est_roll(1) = 0;      % Initializing vehicle attitude
est_pitch(1) = 0;
est_yaw(1) = yaw_prof(1);

laterr=0;  longerr=0;  alphaerr=0;        % INITIALIZATION in this whole section
height = height_prof(1); height_err = 0;
height1 = height_prof(1);  height2 = height_prof(1);
est_lat(1) = lat_prof(1); est_lat(2) = lat_prof(1);
est_lon(1) = lon_prof(1);
vx1 = vel_prof_L(1,1); vx2 = vx1;
vy1 = vel_prof_L(1,2); vy2 = vy1;
vel_l(1,:) = [vx2 vy2 0];
vel2 = [vx1 vy1 0];  vel1 = vel2;
lat2 = lat_prof(1);  lat1 = lat_prof(1);
DCMnb = [DCMnb_prof(1,1:3); DCMnb_prof(1,4:6); DCMnb_prof(1,7:9)];
est_DCMbn = DCMnb';
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

dverr = zeros(size(dvtot));                  % No errors for this demo

est_dv = dvtot + dverr;           % form profile of 'measured' delta-V's

fprintf(1,' . . . . . . . . . \n')
fprintf(1,' Starting nav computations \n')
C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

h = waitbar(0,'Time Loop: Performing INS pos/vel/attitude updating');
for i = 2:npts,
   td12 = time(i) - time(i-1);
   tdex = 0.5*td12;
   tdint = time(i) - time(i-1);
   
   est_DCMbn = bodupdat(est_DCMbn,est_dtheta(i-1,1:3));
   
   [DCM_ll_I, omega_el_L, omega_ie_L] = lclevupd(lat1,lat2,vx1,vx2,vy1,vy2,...
      height1,height2,td12,tdex,tdint,est_DCMel,vertmech,1,earthflg);
   
   est_DCMbn = C*(DCM_ll_I*(C*est_DCMbn));
   
   eul_vect = dcm2eulr(est_DCMbn);
   est_roll(i) = eul_vect(1);
   est_pitch(i) = eul_vect(2);
   est_yaw(i) = eul_vect(3);
   
   est_delv_b = est_dv(i-1,1:3);       % extract delta-V for current point in time

   del_Vl = C*(est_DCMbn*est_delv_b');

   omega1_el_L = omega2_el_L;   omega2_el_L = omega_el_L;
   [est_DCMel, DCM_ll_E] = navupdat(omega1_el_L,omega2_el_L,td12,est_DCMel,1);
   
   h_extr = extrapol(height1,height2,td12,tdex);
   lat_extr = extrapol(lat1,lat2,td12,tdex);
   g_extr = gravity(lat_extr,h_extr);
 
   vtmp = velupdat(vel2,vel1,td12,tdex,del_Vl,...
      omega_el_L,est_DCMel,g_extr,0,tdint);
   
   vel_l(i,:) = vtmp';
   
   est_height(i,1) = height_prof(i);
   height1 = height2;  height2 = est_height(i,1);
    vx1 = vx2; vy1 = vy2;
    vx2 = vel_l(i,1);  vy2 = vel_l(i,2);
    vel1 = vel2;   vel2 = vel_l(i,:);
    llw_vect = dcm2llw(est_DCMel);
    est_lat(i) = llw_vect(1); est_lon(i) = llw_vect(2);  est_alpha(i) = llw_vect(3);
    lat1 = lat2;  lat2 = est_lat(i);
    laterr(i) = est_lat(i) - lat_prof(i);
    longerr(i) = est_lon(i) - lon_prof(i);
    
    waitbar(i/npts,h)
end
close(h)

est_alpha0 = est_alpha;
est_lat0 = est_lat;
est_lon0 = est_lon;
vel_l_0 = vel_l;
est_roll0 = est_roll;
est_pitch0 = est_pitch;
est_yaw0 = est_yaw;

%%save dyn_flt_2_ins_0 est_alpha0 est_lat0 est_lon0 vel_l_0 est_roll0 est_pitch0 est_yaw0 -V6   % Version 7.x
save dyn_flt_2_ins_0 est_alpha0 est_lat0 est_lon0 vel_l_0 est_roll0 est_pitch0 est_yaw0   % Version 6.x

close all
subplot(211)
plot(time/60,lat_prof*180/pi,time/60,est_lat*180/pi)
title('DYN FLIGHT INS - Baseline without sensor or initialization errors')
ylabel('latitude in degrees')
xlabel('time in minutes')

subplot(212)
plot(time/60,lon_prof*180/pi,time/60,est_lon*180/pi)
ylabel('longitude in degrees')
xlabel('time in minutes')

figure
plot(time/60,laterr*(180/pi)*3600,time/60,longerr*(180/pi)*3600)
title('DYN FLIGHT INS:  Algorithm Errors')
ylabel('error in arc-seconds')
xlabel('run time in minutes')
