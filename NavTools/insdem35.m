%
%   insdem35.m
%
%    04/11/05;  08/18/07
%    M. & S. Braasch
%    GPSoft LLC
%

clear all
close all

%insdem32
load dem32_dat

%%%%%%%%%%% Sensor Error Definitions  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%                                    % generate the component of delta-theta
%                                    % associated with craft rate
deltacr = gendelcr2(lat_prof,vel_prof_L,height_prof,...
                       time,DCMnb_prof,DCMel_prof,earthflg);

deltheta = dthetbody + deltaer + deltacr;   % ideal (error-free) gyro output

dtherr = gentherr(deltheta,time,dthparam,98765);  % generate delta-theta errors

est_dtheta = deltheta + dtherr;   % form profiles of 'measured' delta-theta's

roll(1) = 0;      % Initializing vehicle attitude
pitch(1) = 0;
yaw(1) = yaw_prof(1);

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

dverr = gendverr(dvtot,time,dvparam,76543);  % Generate delta-V errors

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
   roll(i) = eul_vect(1);
   pitch(i) = eul_vect(2);
   yaw(i) = eul_vect(3);
   
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

load dem34_0

N = max(size(est_lat0));                           % Compute horizontal position
for i = 1:N,                                       % error by finding the ENU
   truxyz = llh2xyz([est_lat0(i) est_lon0(i) 0]);  % coordinates of the 
   insxyz = llh2xyz([est_lat(i) est_lon(i) 0]);    % INS-derived position 
   enu = xyz2enu(insxyz,truxyz);                   % relative to the truth
   horz_pos_err(i) = norm(enu);
end

close all
subplot(211)
plot(time/60,lat_prof*180/pi,time/60,est_lat*180/pi)
title('INSDEM35')
ylabel('latitude in degrees')
xlabel('time in minutes')

subplot(212)
plot(time/60,lon_prof*180/pi,time/60,est_lon*180/pi)
ylabel('longitude in degrees')
xlabel('time in minutes')
text(0.5,-81.5,'Note that in both plots the true and')
text(0.5,-81.6,'computed routes overlay each other')


subplot(211)
plot(time/60,(est_lat-est_lat0)*180/pi,time/60,(est_lon-est_lon0)*180/pi)
%axis([0 2 -0.06 0.02])
title('INSDEM35 - Position and Velocity Errors')
ylabel('error in degrees')
xlabel('time in minutes')
legend('lat','long','Location','Best')
%text(0.3,-0.02,'LAT')
%text(1.6,-0.03,'LONG')

subplot(212)
plot(time/60,vel_l(:,1:2)-vel_l_0(:,1:2))
ylabel('velocity error in m/s')
xlabel('time in minutes')
%text(0.5,1,'EAST')
%text(0.3,-2,'NORTH')
legend('east','north','Location','Best')

figure
subplot(211)
plot(time/60,(roll-roll0)*180/pi,time/60,(pitch-pitch0)*180/pi)
%axis([0 9 -0.1 0.1])
title('INSDEM35 - Euler Angle Errors')
ylabel('error in degrees')
xlabel('time in minutes')
%text(1.05,-0.01,'ROLL')
%text(0.6,0.01,'PITCH')
legend('roll','pitch','Location','Best')

subplot(212)
plot(time/60,(yaw-yaw0)*180/pi)
ylabel('yaw error in degrees')
xlabel('time in minutes')
%text(1.2,0,'YAW')

figure
%%m2nmi = 3.2808/6076;      % Convert from meters to nautical miles
plot(time/60,horz_pos_err)
title('INSDEM35 - Horizontal Position Error')
ylabel('error in meters')
xlabel('time in minutes')


