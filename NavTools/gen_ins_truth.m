function gen_ins_truth(profilename,earthflg,vertmech,...
    latstart,lonstart,outfilename)
%GEN_INS_TRUTH	Generate 'truth' INS measurements (delta-Vs and
%delta-thetas) and 'truth' INS position, velocity and attitude.
%
%	GEN_INS_TRUTH(profilename,earthflg,vertmech,latstart,lonstart,...
%                 outfilename)	
%
%    INPUTS
%	profilename = filename of MAT-file containing the output of the vehicle
%	              trajectory generator
%   earthflg = earth shape flag
%              0 = spherical earth with radius 6378137 meters
%              1 = ellipsoidal earth (WGS-84)
%   vertmech = Argument specifying vertical mechanization
%              0 = north-pointing (default)
%              1 = free azimuth (vertical craft rate set to zero)
%              2 = foucault (vertical component of spatial rate set to zero)
%              3 = unipolar (wander angle rate set to +/- longitude
%                            rate; - for northern hemisphere;
%                            + for southern hemisphere)
%   latstart = latitude of starting point of vehicle trajectory (deg)
%   lonstart = longitude of starting point of vehicle trajectory (deg)
%   outname = filename given to the MAT-file containing the output of this
%             program
%
%    OUTPUTS
%	NONE.  The results are saved in MAT-file specified by OUTNAME
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
%  September 2009
%  Copyright (c) 2009 by GPSoft LLC
%  All Rights Reserved.
%


load(profilename)    % load the vehicle trajectory file

orgllh = [latstart*pi/180 lonstart*pi/180 0];  % Start-point of ENU profile

npts = length(time);

vel_prof_L = profile(:,4:6);             % extracting velocity profile

DCMnb_prof = profile(:,10:18);           % extracting true nav-to-body
%                                        % direction cosine matrix
fprintf(1,' Generating body-motion component of delta-theta \n')
dthetbody = gendthet(DCMnb_prof);       % component of delta-theta associated
%                                      % with body motion

dtherr = [0 0 0];  % no errors in the baseline run

% Preallocate arrays for speed:
deltheta = zeros(npts-1,3);
dvtot = zeros(npts-1,3);
est_lat = zeros(npts,1);
est_lon = zeros(npts,1);
est_alpha = zeros(npts,1);
est_height = zeros(npts,1);
est_roll = zeros(npts,1);
est_pitch = zeros(npts,1);
est_yaw = zeros(npts,1);
vel_L = zeros(npts,3);

 % INITIALIZATION in this whole section
est_height(1) = orgllh(3);
height1 = orgllh(3);  height2 = orgllh(3);
est_lat(1) = orgllh(1); est_lat(2) = orgllh(1);
est_lon(1) = orgllh(2);
vx1 = vel_prof_L(1,1); vx2 = vx1;
vy1 = vel_prof_L(1,2); vy2 = vy1;
vel_L(1,:) = [vx2 vy2 0];
vel2 = [vx1 vy1 0];  vel1 = vel2;
lat2 = orgllh(1);  lat1 = orgllh(1);
DCMnb = [DCMnb_prof(1,1:3); DCMnb_prof(1,4:6); DCMnb_prof(1,7:9)];
est_DCMbn = DCMnb';
    tmp=dcm2eulr(est_DCMbn);
    est_roll(1) = tmp(1);
    est_pitch(1) = tmp(2);
    est_yaw(1) = tmp(3);

est_DCMel = llw2dcm([orgllh(1) orgllh(2) 0]);

omega2_el_L = crafrate(orgllh(1),vx1,vy1,orgllh(3),est_DCMel,earthflg,vertmech);

fprintf(1,' Generating body-motion component of delta-V \n')
deltav_b = gendv(vel_prof_L,DCMnb_prof);    % Generate the component of delta-V
%                                           % associated with body motion relative
%                                           % to the earth

waitbarincr = floor(npts/10);

fprintf(1,' . . . . . . . . . \n')
fprintf(1,' Starting nav computations \n')
C = [0 1 0; 1 0 0; 0 0 -1];    % conversion between NED and ENU

h = waitbar(0,'Time Loop: Performing INS pos/vel/attitude updating');
for i = 2:npts,
   td12 = time(i) - time(i-1);
   tdex = 0.5*td12;
   tdint = time(i) - time(i-1);
   
   deltaer = earthrot2(tdint,est_DCMel,est_DCMbn);
   deltacr = gendelcr3(est_lat(i-1),vel_L(i-1,:),est_height(i-1),...
       tdint,est_DCMbn,est_DCMel,earthflg);
   deltheta(i-1,1:3) = dthetbody(i-1,1:3) + deltaer + deltacr;   % ideal (error-free) gyro output
   est_dtheta = deltheta(i-1,1:3) + dtherr;   % form profiles of 'measured' delta-theta's

   est_DCMbn = bodupdat(est_DCMbn,est_dtheta(1,1:3));
   
   [DCM_ll_I, omega_el_L, omega_ie_L] = lclevupd(lat1,lat2,vx1,vx2,vy1,vy2,...
      height1,height2,td12,tdex,tdint,est_DCMel,vertmech,1,earthflg);
   
   est_DCMbn = C*(DCM_ll_I*(C*est_DCMbn));
   
   eul_vect = dcm2eulr(est_DCMbn);
   est_roll(i,1) = eul_vect(1);
   est_pitch(i,1) = eul_vect(2);
   est_yaw(i,1) = eul_vect(3);
   
   dvcor = gendvcor3(est_lat(i-1),vel_L(i-1,:),est_height(i-1),...
       tdint,est_DCMbn,est_DCMel,earthflg);  % compute Coriolis component

   dvtot(i-1,1:3) = deltav_b(i-1,1:3) + dvcor;  % ideal error-free delta-v
   est_dv = dvtot(i-1,1:3);   % No errors being simulated in this baseline run
   
   del_Vl = C*(est_DCMbn*est_dv');

   omega1_el_L = omega2_el_L;   omega2_el_L = omega_el_L;
   [est_DCMel, DCM_ll_E] = navupdat(omega1_el_L,omega2_el_L,td12,est_DCMel,1);
   
   h_extr = extrapol(height1,height2,td12,tdex);
   lat_extr = extrapol(lat1,lat2,td12,tdex);
   g_extr = gravity(lat_extr,h_extr);
 
   vtmp = velupdat(vel2,vel1,td12,tdex,del_Vl,...
      omega_el_L,est_DCMel,g_extr,0,tdint);
   
   vel_L(i,:) = vtmp';
   
   est_height(i,1)= est_height(i-1,1) + tdint*mean([vel_L(i,3); vel_L(i-1,3)]);
   height1 = height2;  height2 = est_height(i,1);
   
   vx1 = vx2; vy1 = vy2;
    vx2 = vel_L(i,1);  vy2 = vel_L(i,2);
    vel1 = vel2;   vel2 = vel_L(i,:);
    llw_vect = dcm2llw(est_DCMel);
    est_lat(i,1) = llw_vect(1); est_lon(i,1) = llw_vect(2);  
    est_alpha(i,1) = llw_vect(3);
    lat1 = lat2;  lat2 = est_lat(i);
    
    if rem(i,waitbarincr) == 0,
        waitbar(i/npts,h)
    end
end
close(h)

tru_alpha = est_alpha;
tru_lat = est_lat;
tru_lon = est_lon;
tru_height = est_height;
tru_vel_L = vel_L;
tru_roll = est_roll;
tru_pitch = est_pitch;
tru_yaw = est_yaw;

save(outfilename,'npts','time','tru_alpha','tru_lat','tru_lon',...
    'tru_height','tru_vel_L','tru_roll','tru_pitch','tru_yaw',...
    'deltheta','dvtot','-V6')

