function profile = pro_f16climb(initpos,initvel,heading,duration,deltat)
%PROPITCH       F-16 flight profile sub-generator for
%               climbing/descending flight.  Local-level
%               (i.e., East-North-Up coordinates) version.
%
%  profile = pro_f16climb(initpos,initvel,heading,duration,deltat)
%
%   INPUTS
%       initpos = initial position of vehicle 
%                 (3 ENU cartesian coordinates) (meters)
%       initvel = initial velocity vector (3 ENU cartesian 
%                 coordinates) (m/s)
%       heading = desired heading of this segment (radians)
%       duration = duration of climb in seconds
%       deltat = time increment in seconds
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

vt = norm(initvel)*m2f;
dt = deltat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rollrate=0; pitchrate=0; turnrate=0; coord=0; stab=0;
thtl_init = 0;
el_init = 0;
ail_init = 0;
rdr_init = 0;

x_prev = x;  xd_prev = xd;
%%x = [vt alpha_init beta_init phi theta psi 0 0 0 initpos(2) initpos(1) initpos(3) 0]';
x = [vt 0 0 0 0 0 0 0 0 0 0 0 0]';
% Note climb angle (gamma) is passed in via global variable
[init_cost,final_cost] = trimf16(thtl_init,el_init,ail_init,rdr_init);
thtl_climb = thtl_final;
el_climb = el_final;
ail_climb = ail_final;
rdr_climb = rdr_final;
x_trim_climb = x;
    x(6) = heading;
    [xd] = f16_6dof(0,x,thtl_climb,el_climb,ail_climb,rdr_climb,xcg);
x = x_prev;
%%
x(4:5) = x_trim_climb(4:5);

    t = deltat:deltat:duration; t=t';
    tmp = ones(max(size(t)),1);

    profile(:,1) = initpos(1)*tmp + xd(11)*t*f2m;
    profile(:,2) = initpos(2)*tmp + xd(10)*t*f2m;
    profile(:,3) = initpos(3)*tmp + xd(12)*t*f2m;

    profile(:,4) = xd(11)*f2m*tmp;
    profile(:,5) = xd(10)*f2m*tmp;
    profile(:,6) = xd(12)*f2m*tmp;

    profile(:,7) = NaN*tmp;
    profile(:,8) = NaN*tmp;
    profile(:,9) = NaN*tmp;

    phi = x(4); theta = x(5); psi = x(6);
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
    x(1) = norm(xd(10:12));
    x(4) = phi;
    x(5) = theta;
    x(6) = psi;
    x(13) = x_trim_climb(13);
 
