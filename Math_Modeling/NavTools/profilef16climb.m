function [profile,thtl_end,el_end,ail_end,rdr_end] = ...
    profilef16climb(alt_final,heading,deltat,...
    thtl_start,el_start,ail_start,rdr_start)
%PROFILEF16CLIMB   F-16 flight profile sub-generator for a
%                  climb or descent to a desired altitude.  
%                  Local-level (i.e., East-North-Up coordinates)
%                  version.
%
%  [profile,thtl_end,el_end,ail_end,rdr_end] = ...
%    profilef16climb(alt_final,heading,duration,deltat,...
%    thtl_start,el_start,ail_start,rdr_start)
%
%   INPUTS
%       alt_final = desired altitude in meters
%       heading = desired heading in radians
%       duration = length of flight path segment in seconds
%       deltat = time increment in seconds
%       thtl_start = throttle setting at the beginning of this flight phase
%       el_start = elevator deflection at the beginning of this phase
%       ail_start = aileron deflection at the beginning of this phase
%       rdr_start = rudder deflection at the beginning of this phase
%
%   NOTE: Climb angle is passed in via global variable
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

%	M. & S. Braasch August 2007;  September 2009
%	Copyright (c) 2007-2009 by GPSoft LLC
%	All Rights Reserved.
%

global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha
global thtl_final el_final ail_final rdr_final
global m2f f2m  % meters/feet conversions;  defined in PROFILEF16.M

if nargin<7,error('insufficient number of input arguments'),end

dt = deltat;
vel_fps = x(1);
alt_final_feet = alt_final*m2f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Compute trim conditions for the climb/descent

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
x = [vel_fps 0 0 0 0 0 0 0 0 0 0 x_prev(12) 0]';
[init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
thtl_climb = thtl_final;
el_climb = el_final;
ail_climb = ail_final;
rdr_climb = rdr_final;
x_trim_climb = x;
pitch_climb = x(5);
x = x_prev; xd = xd_prev;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%htrans = waitbar(0,'Generating climb/descent segment');
fprintf(1,' Generating Climb Segment \n')

k = 1; t(k) = 0;
ail = ail_climb;  rdr = rdr_climb; thtl = thtl_climb;
altdiff = x(12) - alt_final_feet;
altdiff_accum = altdiff;
altdiff_init = abs(altdiff);
Kp_pitch = -0.1*pi/180;
Ki_pitch = -0.01*pi/180;
%Kd_pitch = 
pitchdiff = x_trim_climb(5) - x(5);
pitchdiff_init = pitchdiff;
pitchdiff_accum = pitchdiff*dt;
Kp_el = -10;
Ki_el = -1;
Kd_el = 0;
rolldiff = x_trim_climb(4) - x(4);
rolldiff_accum = rolldiff*dt;
Kp_ail = -100;
Ki_ail = -1;
Kd_ail = 25;

Kp_rdr = 1;

speeddiff = x_trim_climb(1) - x(1);
speeddiff_accum = 0;
Kp_thtl = 1e-4;
Ki_thtl = 1e-4;
thtl_cnt = 0;

while abs(altdiff) > 20, 
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
            thtl = thtl_climb;
        %
        rdr = Kp_rdr*x(9);
        if rdr > 30, rdr = 30; end    
        if rdr < -30, rdr = -30; end
        %
    
    [x] = rk4_f16(t(k),dt,x,thtl,el,ail,rdr,xcg);

    altdiff = x(12) - alt_final*m2f;
    altdiff_accum = altdiff_accum + altdiff*dt;
        pitch_desired = pitch_climb;
    pitchdiff = pitch_desired - x(5);
    pitchdiff_accum = pitchdiff_accum + pitchdiff*dt;
    rolldiff = x_trim_level(4) - x(4);
    rolldiff_accum = rolldiff_accum + rolldiff*dt;

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
     
    %waitbar((1-abs(altdiff)/altdiff_init),htrans)
end
%close(htrans)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

thtl_end = thtl;
el_end = el;
ail_end = ail;
rdr_end = rdr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


