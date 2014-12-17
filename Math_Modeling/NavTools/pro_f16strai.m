function profile = pro_f16strai(initpos,initvel,heading,duration,v_final,deltat)
%PROSTRAI       F-16 flight profile sub-generator for a straight-and-level,
%               constant acceleration or constant velocity flight segment.
%       
%	profile = pro_f16strai(initpos,initvel,heading,duration,v_final,deltat)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (3 ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (3 ENU cartesian 
%                 coordinates) (m/s)
%       heading = desired heading angle on this segment (radians)
%       duration = length of flight segment (seconds)
%       v_final = total final velocity (i.e., speed) in meters-per-second
%       deltat = time increment in seconds
%
%   OUTPUTS
%       profile = flight profile
%          profile(i,1:3) = ENU path generated; 1=x, 2=y, 3=z;
%                           units are meters
%          profile(i,4:6) = ENU velocity; 4 = x-velocity,
%                           5 = y-velocity, 6 = z-velocity;
%                           units are meters/second
%          profile(i,7:9) = NOT USED (format retained for compatibility
%                           with other toolbox functions
%          profile(i,10:18) = elements of the direction cosine matrix
%                            (DCM) for vehicle attitude; 10 = DCM(1,1),
%                            11 = DCM(1,2), 12 = DCM(1,3),
%                            13 = DCM(2,1), et cetera
%

%	M. & S. Braasch 02-2005
%	Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%
global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha
global thtl_final el_final ail_final rdr_final

if nargin<6,error('insufficient number of input arguments'),end

f2m = 0.3048;  m2f = 1/f2m;   % feet/meters conversions
% dcmbn=initdcm';
% eulvect=dcm2eulr(dcmbn);
% phi=eulvect(1);
% theta=eulvect(2);
% psi=eulvect(3);
vt_init = norm(initvel)*3.2808;
dt = deltat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gamma = 0;
rollrate=0; pitchrate=0; turnrate=0; coord=0; stab=0;
thtl_init = 0;
el_init = 0;
ail_init = 0;
rdr_init = 0;
alpha_init = 0;
beta_init = 0;

if vt_init >= 300,
    x = [vt_init 0 0 0 0 0 0 0 0 0 0 0 0]';
    [init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
    theta = x(5);
else
    theta = 0;
end

x_prev = x; xd_prev = xd;
vt_fps=v_final*m2f;
%%x = [vt_fps alpha_init beta_init phi theta psi 0 0 0 initpos(2) initpos(1) initpos(3) 0]';
x = [vt_fps 0 0 0 0 0 0 0 0 0 0 0 0]';
[init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
thtl_strai = thtl_final;
el_strai = el_final;
ail_strai = ail_final;
rdr_strai = rdr_final;
x_trim_strai = x;
x = x_prev; xd = xd_prev;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = deltat:deltat:duration; t=t';
    tmp = ones(max(size(t)),1);
    
    %dcmbn = initdcm';
    %eulr_vec = dcm2eulr(dcmbn);
    gnd_trk_ang = heading;
    v_init = norm(initvel);
    accel = (v_final - v_init)/duration;

    profile(:,7) = NaN*tmp;
    profile(:,8) = NaN*tmp;
    profile(:,9) = NaN*tmp;

%     profile(:,4) = initvel(1)*tmp + accel*sin(gnd_trk_ang)*t;
%     profile(:,5) = initvel(2)*tmp + accel*cos(gnd_trk_ang)*t;
%     profile(:,6) = 0;
% 
%     profile(:,1) = initpos(1)*tmp + initvel(1)*t + 0.5*accel*sin(gnd_trk_ang)*t.*t;
%     profile(:,2) = initpos(2)*tmp + initvel(2)*t + 0.5*accel*cos(gnd_trk_ang)*t.*t;
%     profile(:,3) = initpos(3)*tmp;

    profile(:,4) = v_init*sin(gnd_trk_ang)*tmp + accel*sin(gnd_trk_ang)*t;
    profile(:,5) = v_init*cos(gnd_trk_ang)*tmp + accel*cos(gnd_trk_ang)*t;
    profile(:,6) = 0;

    profile(:,1) = initpos(1)*tmp + v_init*sin(gnd_trk_ang)*t + 0.5*accel*sin(gnd_trk_ang)*t.*t;
    profile(:,2) = initpos(2)*tmp + v_init*cos(gnd_trk_ang)*t + 0.5*accel*cos(gnd_trk_ang)*t.*t;
    profile(:,3) = initpos(3)*tmp;

    phi = 0;
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
    x(13) = x_trim_strai(13);
    