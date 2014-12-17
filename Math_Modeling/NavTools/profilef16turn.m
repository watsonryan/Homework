function [profile,thtl_end,el_end,ail_end,rdr_end] = ...
          profilef16turn(hdg_final,altitude,deltat,thtl_start,el_start,...
          ail_start,rdr_start);
%PROFILEF16TURN    F-16 flight profile sub-generator for a coordinated turn.  
%               Local-level (i.e., East-North-Up coordinates) version.
%       
%  [profile,thtl_end,el_end,ail_end,rdr_end] = ...
%          profilef16turn(hdg_final,altitude,deltat,thtl_start,el_start,...
%          ail_start,rdr_start);
%
%   INPUTS
%       hdg_final = desired final heading angle in radians
%       altitude = desired altitude in meters
%       deltat = time increment in seconds
%       thtl_start = throttle setting at the beginning of this flight phase
%       el_start = elevator deflection at the beginning of this phase
%       ail_start = aileron deflection at the beginning of this phase
%       rdr_start = rudder deflection at the beginning of this phase
%
%   NOTE: turn rate is passed in via global variable
%
%   OUTPUTS
%       profile = flight profile
%          profile(i,1:3) = ENU path generated; 1=x, 2=y, 3=z
%          profile(i,4:6) = ENU velocity; 4 = x-velocity,
%                           5 = y-velocity, 6 = z-velocity 
%          profile(i,7:9) = NOT USED
%          profile(i,10:18) = elements of the direction cosine matrix
%                            (DCM) for vehicle attitude; 10 = DCM(1,1),
%                            11 = DCM(1,2), 12 = DCM(1,3),
%                            13 = DCM(2,1), et cetera
%       thtl_end = throttle setting at the end of this flight phase
%       el_end = elevator deflection at the end of this phase
%       ail_end = aileron deflection at the end of this phase
%       rdr_end = rudder deflection at the end of this phase
%

%   Waitbars commented out to increase speed of execution (Sept 2009)

%	M. & S. Braasch 02-05;  08-2007; 09-2009
%	Copyright (c) 2005-2009 by GPSoft LLC
%	All Rights Reserved.
%

global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha
global thtl_final el_final ail_final rdr_final

if nargin<7,error('insufficient number of input arguments'),end

f2m = 0.3048;  m2f = 1/f2m;   % feet/meters conversions

dt = deltat;
vel_fps = x(1);
alt_final_feet = altitude*m2f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Compute trim conditions for the turn

    centrip_acc = vel_fps*turnrate/32.2;
    phi_init = atan(centrip_acc);
    Q_init = sin(phi_init)*turnrate;   % equation 1.3-21
    R_init = cos(phi_init)*turnrate;
gamma=0;
rollrate=0; pitchrate=0; coord=1; stab=0;
thtl_init = 0;
el_init = 0;
ail_init = 0;
rdr_init = 0;
alpha_init = 0;
beta_init = 0;

x_prev = x;
xd_prev = xd;
% gamma is passed in via global variable
x = [vel_fps 0 0 0 0 0 0 Q_init R_init 0 0 alt_final_feet 0]';
[init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
thtl_turn = thtl_final;
el_turn = el_final;
ail_turn = ail_final;
rdr_turn = rdr_final;
x_trim_turn = x;
roll_turn = x(4);
pitch_turn = x(5);
pitch_rate = x(8);
x = x_prev; xd = xd_prev;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute trim conditions for level flight
gamma_old = gamma;
gamma = 0;

rollrate=0; pitchrate=0; turnrate=0; coord=0; stab=0;
thtl_init = 0;
el_init = 0;
ail_init = 0;
rdr_init = 0;
alpha_init = 0;
beta_init = 0;

x_prev = x;
xd_prev = xd;
% gamma is passed in via global variable
x = [vel_fps 0 0 0 0 0 0 0 0 0 0 alt_final_feet 0]';
[init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
thtl_level = thtl_final;
el_level = el_final;
ail_level = ail_final;
rdr_level = rdr_final;
x_trim_level = x;
pitch_level = x(5);
x = x_prev; xd = xd_prev;

gamma = gamma_old;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%htrans = waitbar(0,'Generating turn segment');
fprintf(1,' Generating Turn Segment \n')

k = 1; t(k) = 0;
ail = ail_level;  rdr = rdr_level; thtl = thtl_level; el = el_level;
altdiff = 0;
altdiff_accum = altdiff;
altdiff_init = abs(altdiff);
Kp_pitch = -0.1*pi/180;
Ki_pitch = -0.01*pi/180;
%Kd_pitch = 
pitchdiff = x_trim_turn(5) - x(5);
pitchdiff_init = pitchdiff;
pitchdiff_accum = pitchdiff*dt;
Kp_el = -200;
Ki_el = -20;
Kd_el = 0;
rolldiff = x_trim_turn(4) - x(4);
rolldiff_accum = rolldiff*dt;
Kp_ail = -100;
Ki_ail = -1;
Kd_ail = 25;

hdg_curr = x(6);
  if hdg_curr > 2*pi, hdg_curr = hdg_curr - 2*pi; end
  if hdg_curr < 0, hdg_curr = hdg_curr + 2*pi; end
hdg_diff = hdg_final - hdg_curr;
hdg_diff_init = hdg_diff;
Kp_rdr = 1;

speeddiff = x_trim_turn(1) - x(1);
speeddiff_accum = 0;
Kp_thtl = 1e-4;
Ki_thtl = 1e-4;
thtl_cnt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = x(6);
if psi < 0, psi = psi + 2*pi; end
if psi > 2*pi, psi = psi - 2*pi; end
turndelta = 1*pi/180;   
psi_final = hdg_final;
if psi_final>2*pi, psi_final=psi_final-2*pi; end
ang_diff = psi_final - psi;
if (turnrate<0)&&(ang_diff>0), ang_diff=ang_diff-2*pi; end
if (turnrate>0)&&(ang_diff<0), ang_diff=ang_diff+2*pi; end
ang_diff_init = abs(ang_diff);
prev_mag_ang_diff = ang_diff_init;
flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while flag == 0,
    k = k + 1;
    t(k) = t(k-1) + dt;
        el = Kp_el*pitchdiff + Ki_el*pitchdiff_accum + Kd_el*x(8);
        if el > 25, el = 25; end
        if el < -25, el = -25; end
        %
        ail = Kp_ail*rolldiff + Ki_ail*rolldiff_accum + Kd_ail*x(7);
            if ail > 21.5, ail = 21.5; end
            if ail < -21.5, ail = -21.5; end
        %
            thtl = thtl_turn;
        %
        rdr = Kp_rdr*(x(9)-turnrate);   
        if rdr > 30, rdr = 30; end    
        if rdr < -30, rdr = -30; end
        %
    
    [x] = rk4_f16(t(k),dt,x,thtl,el,ail,rdr,xcg);

    altdiff = x(12) - alt_final_feet;
    altdiff_accum = altdiff_accum + altdiff*dt;
    alt_pitch_desired = Kp_pitch*altdiff + Ki_pitch*altdiff_accum;
    pitchdiff = x_trim_turn(5) - x(5) + alt_pitch_desired;
    pitchdiff_accum = pitchdiff_accum + pitchdiff*dt;
    rolldiff = x_trim_turn(4) - x(4);
    rolldiff_accum = rolldiff_accum + rolldiff*dt;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    psi = x(6);
    if psi < 0, psi = psi + 2*pi; end
    if psi > 2*pi, psi = psi - 2*pi; end

    ang_diff = psi_final - psi;
    if (turnrate<0)&&(ang_diff>0), ang_diff=ang_diff-2*pi; end
    if (turnrate>0)&&(ang_diff<0), ang_diff=ang_diff+2*pi; end

    if abs(ang_diff) < turndelta, flag = 1; end
    prev_mag_ang_diff = abs(ang_diff);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    [xd] = f16_6dof(t(k),x,thtl,el,ail,rdr,xcg);

    i = k - 1;
    
    profile(i,1) = x(11)*f2m; 
    profile(i,2) = x(10)*f2m; 
    profile(i,3) = x(12)*f2m;
    profile(i,4) = xd(11)*f2m; 
    profile(i,5) = xd(10)*f2m; 
    profile(i,6) = xd(12)*f2m;
    profile(i,7:9) = [NaN NaN NaN];
    
    phi = x(4); theta = x(5); psi = x(6);
    DCMnb = eulr2dcm([phi theta psi]);
    dcm = DCMnb;
    profile(i,10) = dcm(1,1);
    profile(i,11) = dcm(1,2);
    profile(i,12) = dcm(1,3);
    profile(i,13) = dcm(2,1);
    profile(i,14) = dcm(2,2);
    profile(i,15) = dcm(2,3);
    profile(i,16) = dcm(3,1);
    profile(i,17) = dcm(3,2);
    profile(i,18) = dcm(3,3);
    
    %waitbar(1-abs(ang_diff)/ang_diff_init,htrans)

end
%close(htrans)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%htrans = waitbar(0,'Generating transition to level segment');
fprintf(1,' Generating Transition to Level Segment \n')

altdiff = 0;
altdiff_accum = altdiff;
altdiff_init = abs(altdiff);
Kp_pitch = -0.1*pi/180;
Ki_pitch = -0.01*pi/180;
%Kd_pitch = 
pitchdiff = x_trim_turn(5) - x(5);
pitchdiff_init = pitchdiff;
pitchdiff_accum = pitchdiff*dt;
Kp_el = -200;
Ki_el = -20;
Kd_el = 0;
rolldiff = x_trim_turn(4) - x(4);
rolldiff_accum = rolldiff*dt;
Kp_ail = -100;
Ki_ail = -1;
Kd_ail = 25;

hdg_curr = x(6);
  if hdg_curr > 2*pi, hdg_curr = hdg_curr - 2*pi; end
  if hdg_curr < 0, hdg_curr = hdg_curr + 2*pi; end
hdg_diff = hdg_final - hdg_curr;
hdg_diff_init = hdg_diff;
Kp_rdr = 1;

speeddiff = x_trim_turn(1) - x(1);
speeddiff_accum = 0;
Kp_thtl = 1e-4;
Ki_thtl = 1e-4;
thtl_cnt = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = x(6);
if psi < 0, psi = psi + 2*pi; end
if psi > 2*pi, psi = psi - 2*pi; end
turndelta = 0.1*pi/180;   
psi_final = hdg_final;
if psi_final>2*pi, psi_final=psi_final-2*pi; end
ang_diff = psi_final - psi;
if (turnrate<0)&&(ang_diff>0), ang_diff=ang_diff-2*pi; end
if (turnrate>0)&&(ang_diff<0), ang_diff=ang_diff+2*pi; end
ang_diff_init = abs(ang_diff);
prev_mag_ang_diff = ang_diff_init;
flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while flag == 0,
    k = k + 1;
    t(k) = t(k-1) + dt;
        el = Kp_el*pitchdiff + Ki_el*pitchdiff_accum + Kd_el*x(8);
        if el > 25, el = 25; end
        if el < -25, el = -25; end
        %
        ail = Kp_ail*rolldiff + Ki_ail*rolldiff_accum + Kd_ail*x(7);
            if ail > 21.5, ail = 21.5; end
            if ail < -21.5, ail = -21.5; end
        %
            thtl = thtl_turn;
        %
        rdr = Kp_rdr*x(9) - 100*ang_diff;   % use rudder only to make
                                            % final heading adjustment
        if rdr > 30, rdr = 30; end     
        if rdr < -30, rdr = -30; end
        %
    
    [x] = rk4_f16(t(k),dt,x,thtl,el,ail,rdr,xcg);

    altdiff = x(12) - alt_final_feet;
    altdiff_accum = altdiff_accum + altdiff*dt;
    alt_pitch_desired = Kp_pitch*altdiff + Ki_pitch*altdiff_accum;
    pitchdiff = x_trim_level(5) - x(5) + alt_pitch_desired;
    pitchdiff_accum = pitchdiff_accum + pitchdiff*dt;
    rolldiff = x_trim_level(4) - x(4);
    rolldiff_accum = rolldiff_accum + rolldiff*dt;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    psi = x(6);
    if psi < 0, psi = psi + 2*pi; end
    if psi > 2*pi, psi = psi - 2*pi; end

    ang_diff = psi_final - psi;
    if (turnrate<0)&&(ang_diff>0), ang_diff=ang_diff-2*pi; end
    if (turnrate>0)&&(ang_diff<0), ang_diff=ang_diff+2*pi; end

    if abs(ang_diff) < turndelta, flag = 1; end
    prev_mag_ang_diff = abs(ang_diff);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    [xd] = f16_6dof(t(k),x,thtl,el,ail,rdr,xcg);

    i = k - 1;
    
    thtl_save(i) = thtl;
    el_save(i) = el;
    ail_save(i) = ail;
    rdr_save(i) = rdr;
    
    profile(i,1) = x(11)*f2m; 
    profile(i,2) = x(10)*f2m; 
    profile(i,3) = x(12)*f2m;
    profile(i,4) = xd(11)*f2m; 
    profile(i,5) = xd(10)*f2m; 
    profile(i,6) = xd(12)*f2m;
    profile(i,7:9) = [NaN NaN NaN];
    
    phi = x(4); theta = x(5); psi = x(6);
    DCMnb = eulr2dcm([phi theta psi]);
    dcm = DCMnb;
    profile(i,10) = dcm(1,1);
    profile(i,11) = dcm(1,2);
    profile(i,12) = dcm(1,3);
    profile(i,13) = dcm(2,1);
    profile(i,14) = dcm(2,2);
    profile(i,15) = dcm(2,3);
    profile(i,16) = dcm(3,1);
    profile(i,17) = dcm(3,2);
    profile(i,18) = dcm(3,3);
    
    %waitbar(1-abs(ang_diff)/ang_diff_init,htrans)

end
%close(htrans)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

thtl_end = thtl;
el_end = el;
ail_end = ail;
rdr_end = rdr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


