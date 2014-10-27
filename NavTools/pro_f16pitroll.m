function [profile,errflg] = ...
    pro_f16pitroll(initpos,initvel,initdcm,deltat)
%PROPITCH       F-16 flight profile sub-generator for a transition
%               between level flight and a coordinated turn.  
%               Local-level (i.e., East-
%               North-Up coordinates) version.
%       
%  [profile,errflg] = ...
%          pro_f16pitroll(initpos,initvel,initdcm,deltat)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (3 ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (3 ENU cartesian 
%                 coordinates) (m/s)
%       initdcm = initial direction cosine matrix for
%                 vehicle attitude (navigation to body
%                 frame) (3x3 matrix)
%       deltat = time increment in seconds
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
%        errflg = 0 if there is no error condition
%               = 1 if pitch manuever cannot be simulated and the user should
%                 specify a smaller time step or smaller pitch rate
%

%	M. & S. Braasch 02-05
%	Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha
global thtl_final el_final ail_final rdr_final

if nargin<4,error('insufficient number of input arguments'),end
[m,n]=size(initpos); if m>n, pos=initpos'; else, pos=initpos; end
[m,n]=size(initvel); if m>n, vel=initvel'; else, vel=initvel; end

f2m = 0.3048;  m2f = 1/f2m;   % feet/meters conversions

errflg = 0;
dcmbn=initdcm';
eulvect=dcm2eulr(dcmbn);
phi=eulvect(1); theta=eulvect(2); psi=eulvect(3);
vt = norm(initvel)*m2f;
dt = deltat;
prevpos = pos;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    centrip_acc = vt*turnrate/32.2;
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
%%x = [vt alpha_init beta_init phi theta psi 0 0 0 initpos(2) initpos(1) initpos(3) 0]';
x = [vt 0 0 0 0 0 0 Q_init R_init 0 0 0 0]';
[init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
thtl_turn = thtl_final;
el_turn = el_final;
ail_turn = ail_final;
rdr_turn = rdr_final;
x_trim_turn = x;
x = x_prev; xd = xd_prev;
%%x(10:12) = x_prev(10:12); xd(10:12) = xd_prev(10:12);
%%x(6) = x_prev(6); xd(6) = xd_prev(6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k = 1; t(k) = 0;
thtl = thtl_turn;
rdr = rdr_turn;
%htrans = waitbar(0,'Generating attitude transition segment');
%trim_flag = 0;  coord = 0; rollrate=pi; pitchrate=0; turnrate=0;
pitchdiff = x_trim_turn(5) - x(5);
rolldiff = x_trim_turn(4) - x(4);
srdiff = sign(rolldiff);
while sign(rolldiff) == srdiff, %NOTE!! ASSUMES PITCH CONVERGES BEFORE ROLL
    k = k + 1;
    t(k) = t(k-1) + dt;
        %[init_cost,final_cost] = trim_rt2_f16(thtl,el,ail,rdr);
        %thtl = thtl_turn;
        el = el_turn - 100*pitchdiff;
            if el > 25, el = 25; end
            if el < -25, el = -25; end
        ail = -20*rolldiff/(pi/2);   % negative ail produces positive roll
            if ail > 21.5, ail = 21.5; end
            if ail < -21.5, ail = -21.5; end
        %rdr = rdr_turn;
        %
        %
    [x] = rk4_f16(t(k),dt,x,thtl,el,ail,rdr,xcg);
    pitchdiff = x_trim_turn(5) - x(5);
    rolldiff = x_trim_turn(4) - x(4);

    [xd] = f16_6dof(t(k),x,thtl,el,ail,rdr,xcg);

    i = k - 1;
    
    %profile(i,1:3) = x(10:12)';
    profile(i,1) = x(11)*f2m; 
    profile(i,2) = x(10)*f2m; 
    profile(i,3) = x(12)*f2m;
    %profile(i,4:6) = xd(10:12)';
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
    
    %waitbar(x(4)/x_trim_turn(4),htrans)
end
%close(htrans)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
