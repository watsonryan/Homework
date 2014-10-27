function [profile,thtl_end,el_end,ail_end,rdr_end] = ...
          profilef16start(v_final,deltat,duration)
%PROFILEF16START   F-16 flight profile sub-generator for initial 
%                  acceleration on the ground (runway).  Local-level 
%                  (i.e., East-North-Up coordinates) version.
%       
%       [profile,thtl_end,el_end,ail_end,rdr_end] = ...
%           profilef16start(v_final,deltat,duration);
%
%   INPUTS
%       v_final = desired airspeed in meters-per-second
%       deltat = time increment in seconds
%       duration = length of the ground acceleration segment in seconds
%
%   NOTE: This function is used to simulate the acceleration of the
%         vehicle on the runway prior to rotation.  Pitch and roll
%         are thus kept at zero.
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

%	M. & S. Braasch September 2007
%	Copyright (c) 2007 by GPSoft LLC
%	All Rights Reserved.
%

global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha
global thtl_final el_final ail_final rdr_final
global m2f f2m   % meters/feet conversion; defined in PROFILEF16.M

if nargin<3,error('insufficient number of input arguments'),end

dt = deltat;

psi_constant = x(6);
v_final_fps = v_final*m2f;
v_init = x(1)*f2m;

initpos(1) = x(11)*f2m;
initpos(2) = x(10)*f2m;
initpos(3) = x(12)*f2m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 0;
rollrate=0; pitchrate=0; turnrate=0; coord=0; stab=0;
thtl_init = 0;
el_init = 0;
ail_init = 0;
rdr_init = 0;
alpha_init = 0;
beta_init = 0;

x_prev = x; xd_prev = xd;
x = [v_final_fps 0 0 0 0 0 0 0 0 0 0 0 80]';
[init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
thtl_start = thtl_final;
el_start = el_final;
ail_start = ail_final;
rdr_start = rdr_final;
x_trim_final = x;
x = x_prev; xd = xd_prev;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = deltat:deltat:duration; t=t';
    tmp = ones(max(size(t)),1);
    
    gnd_trk_ang = psi_constant;
    accel = (v_final - v_init)/duration;

    profile(:,7) = NaN*tmp;
    profile(:,8) = NaN*tmp;
    profile(:,9) = NaN*tmp;

    profile(:,4) = v_init*sin(gnd_trk_ang)*tmp + accel*sin(gnd_trk_ang)*t;
    profile(:,5) = v_init*cos(gnd_trk_ang)*tmp + accel*cos(gnd_trk_ang)*t;
    profile(:,6) = 0;

    profile(:,1) = initpos(1)*tmp + v_init*sin(gnd_trk_ang)*t + 0.5*accel*sin(gnd_trk_ang)*t.*t;
    profile(:,2) = initpos(2)*tmp + v_init*cos(gnd_trk_ang)*t + 0.5*accel*cos(gnd_trk_ang)*t.*t;
    profile(:,3) = initpos(3)*tmp;

    phi = 0;
    theta = 0;
    psi = gnd_trk_ang;
    DCMnb = eulr2dcm([phi theta psi]);
    dcm = DCMnb;
    profile(:,10) = dcm(1,1)*tmp;
    profile(:,11) = dcm(1,2)*tmp;
    profile(:,12) = dcm(1,3)*tmp;
    profile(:,13) = dcm(2,1)*tmp;
    profile(:,14) = dcm(2,2)*tmp;
    profile(:,15) = dcm(2,3)*tmp;
    profile(:,16) = dcm(3,1)*tmp;
    profile(:,17) = dcm(3,2)*tmp;
    profile(:,18) = dcm(3,3)*tmp;

    k = size(profile,1);
    x(10) = profile(k,2)*m2f;
    x(11) = profile(k,1)*m2f;
    x(12) = profile(k,3)*m2f;
    xd(10) = profile(k,5)*m2f;
    xd(11) = profile(k,4)*m2f;
    xd(12) = profile(k,6)*m2f;
    x(1) = norm(xd(10:12));
    x(4) = phi;
    x(5) = theta;
    x(6) = psi;
    x(13) = x_trim_final(13);
    

thtl_end = thtl_start;
el_end = el_start;
ail_end = ail_start;
rdr_end = rdr_start;
