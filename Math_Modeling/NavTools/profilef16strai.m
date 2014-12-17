function [profile,thtl_end,el_end,ail_end,rdr_end] = ...
          profilef16strai(altitude,heading,duration,v_final,deltat,...
          thtl_start,el_start,ail_start,rdr_start);
%PROFILEF16STRAI   F-16 flight profile sub-generator for a
%                  level flight segment (i.e., altitude-hold mode).  
%                  Local-level (i.e., East-North-Up coordinates) version.
%       
%       [profile,thtl_end,el_end,ail_end,rdr_end] = ...
%           profilef16strai(altitude,heading,duration,v_final,deltat,...
%           thtl_start,el_start,ail_start,rdr_start);
%
%   INPUTS
%       altitude = desired altitude in meters
%       heading = desired heading in radians
%       duration = length of flight path segment in seconds
%       v_final = desired airspeed in meters-per-second
%       deltat = time increment in seconds
%       thtl_start = throttle setting at the beginning of this flight phase
%       el_start = elevator deflection at the beginning of this phase
%       ail_start = aileron deflection at the beginning of this phase
%       rdr_start = rudder deflection at the beginning of this phase
%
%   NOTE: Climb angle (zero for this segment) is passed in via 
%         global variable
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

%	M. & S. Braasch January 2007; September 2009
%	Copyright (c) 2007-2009 by GPSoft LLC
%	All Rights Reserved.
%

global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha
global thtl_final el_final ail_final rdr_final
global m2f f2m   % meters/feet conversion; defined in PROFILEF16.M

if nargin<9,error('insufficient number of input arguments'),end
if isnan(duration),
    duration_flag = 0;
else
    duration_flag = 1;
end

dt = deltat;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Compute trim conditions for level flight with v_final
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
x = [v_final*m2f 0 0 0 0 0 0 0 0 0 0 altitude*m2f 0]';
[init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
thtl_level = thtl_final;
el_level = el_final;
ail_level = ail_final;
rdr_level = rdr_final;
x_trim_level = x;
x = x_prev; xd = xd_prev;

gamma = gamma_old;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%htrans = waitbar(0,'Generating straight segment');
fprintf(1,' Generating Straight Segment \n')

k = 1; t(k) = 0;
ail = ail_level;  rdr = rdr_level; thtl = thtl_level; el = el_level;
altdiff = x(12) - altitude*m2f;
altdiff_accum = altdiff*dt;
Kp_pitch = -0.1*pi/180;
Ki_pitch = -0.01*pi/180;
%Kd_pitch = 
rolldiff = x_trim_level(4) - x(4);
rolldiff_accum = rolldiff*dt;
Kp_ail = -100;
Ki_ail = -1;
Kd_ail = 25;
pitchdiff = x_trim_level(5) - x(5);
pitchdiff_accum = pitchdiff*dt;
Kp_el = -200;
Ki_el = -20;
Kd_el = 25;
speeddiff = x_trim_level(1) - x(1);
speeddiff_accum = 0;  %speeddiff*dt;
Kp_thtl = 1e-4;
Ki_thtl = 1e-4;
thtl_cnt = 0;
%
Kp_rdr = 1;
%
continueloop = 1;

while continueloop,
    k = k + 1;
    t(k) = t(k-1) + dt;
        %
        el = Kp_el*pitchdiff + Ki_el*pitchdiff_accum + Kd_el*x(8);
        if el > 25, el = 25; end
        if el < -25, el = -25; end
        %
        ail = Kp_ail*rolldiff + Ki_ail*rolldiff_accum + Kd_ail*x(7);
        if ail > 21.5, ail = 21.5; end
        if ail < -21.5, ail = -21.5; end
        %
        % NOTE: In the the throttle conrol section below we are trying
        % to deal with the inherent lag between a change in the throttle
        % position and the resulting steady-state change in speed.  Even
        % though this is a simplified F-16 model, it does a decent job of
        % modeling the lag in the engine response and, more importantly,
        % the finite amount of time it takes for the airplane to move
        % faster or slower in response to the change in engine thrust.
        % Since this lag is considerable (i.e., it may take 10 or 20
        % seconds to reach the new steady-state speed), we do a couple of
        % things.  First, we hold a constant throttle (at fairly high or
        % low level) until the speed gets close to the desired value.
        % Second, after we get close, we only adjust the throttle once
        % every 5 seconds.
        speeddiff = x_trim_level(1) - x(1);
        if speeddiff > 25,
            thtl = 0.8;
        elseif speeddiff < -25,
            thtl = 0.1*thtl_level;
        else
            thtl_cnt = thtl_cnt + dt;
            if thtl_cnt >= 5,
                speeddiff = x_trim_level(1) - x(1);
                speeddiff_accum = speeddiff_accum + speeddiff;  
                thtl = thtl_level + Kp_thtl*speeddiff + Ki_thtl*speeddiff_accum;
                %thtl_save(k) = thtl;
                if thtl < 0, thtl = 0; end
                if thtl > 1, thtl = 1; end
                thtl_cnt = 0;
            end
        end  % end "if speeddiff > 25"
        %
        psi = x(6);
        if psi < 0, psi = psi + 2*pi; end
        if psi > 2*pi, psi = psi - 2*pi; end
        ang_diff = heading - psi;
        if ang_diff > pi, ang_diff = ang_diff - 2*pi; end
        if ang_diff < -pi, ang_diff = ang_diff + 2*pi; end
        rdr = Kp_rdr*x(9) - 100*ang_diff;   % use rudder only to make
                                            % final heading adjustment
        if rdr > 30, rdr = 30; end   
        if rdr < -30, rdr = -30; end
        %
    
    [x] = rk4_f16(t(k),dt,x,thtl,el,ail,rdr,xcg);
    altdiff = x(12) - altitude*m2f;
    altdiff_accum = altdiff_accum + altdiff*dt;
    pitch_desired = Kp_pitch*altdiff + Ki_pitch*altdiff_accum;
    pitchdiff = pitch_desired - x(5);
    pitchdiff_accum = pitchdiff_accum + pitchdiff*dt;
    rolldiff = x_trim_level(4) - x(4);
    rolldiff_accum = rolldiff_accum + rolldiff*dt;

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
    
    if duration_flag,
        continueloop = ( t(k) <= duration );
        %waitbar((t(k)/duration),htrans)
    else
        continueloop = ( (x(1)+5) < v_final*m2f );  
        %waitbar(x(1)/(v_final*m2f),htrans)
    end
    
end
%close(htrans)

thtl_end = thtl;
el_end = el;
ail_end = ail;
rdr_end = rdr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
