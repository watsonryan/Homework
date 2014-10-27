function [xd,an,alat,qbar,amach,q,alpha] = ...
         f16_6dof(time,x,thtl,el,ail,rdr,xcg)
%F16_6DOF     Function to emulate simplified six-degrees-of-freedom
%             non-linear model for an F-16
%
% [xd,an,alat,qbar,amach,q,alpha] = f16_6dof(time,x,thtl,el,ail,rdr,xcg)
%
%  INPUTS
%    time = run time in seconds
%    x is the state vector
%      x(1) = airspeed in feet/sec
%      x(2) = alpha (angle-of-attack) in radians
%      x(3) = beta (side-slip) in radians
%      x(4) = phi (roll angle) in radians
%      x(5) = theta (pitch angle) in radians
%      x(6) = psi (yaw angle) in radians
%      x(7) = P  roll-rate component of aircraft angular velocity vector
%      x(8) = Q  pitch-rate component of aircraft angular velocity vector
%      x(9) = R  yaw-rate component of aircraft angular velocity vector
%      x(10) = north displacement in feet
%      x(11) = east displacement in feet
%      x(12) = altitude in feet
%      x(13) = power setting in percent
%    thtl = throttle setting ranging from 0 to 1
%    el = elevator deflection in degrees
%    ail = aileron deflection in degrees
%    rdr = rudder deflection in degrees
%    xcg = center of gravity position relative to the leading edge of the
%          wing in terms of the mean aerodynamic chord (mac); nominal is
%          0.35mac; a forward cg is 0.3mac and an aft cg is 0.38mac
%
%  PRIMARY OUTPUT
%    xd = derivative of the state vector
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

xd = zeros(13,1);

axx = 9496;
ayy = 55814;
azz = 63100;
axz = 982;  axzs = axz^2;
xpq = axz*(axx-ayy+azz);
gam = axx*azz-axzs;
xqr = azz*(azz-ayy)+axzs;
zpq = (axx-ayy)*axx+axzs;
ypr = azz - axx;
%weight = 25000;
gd = 32.17;
%%mass = weight/gd;
mass = 636.94;

s = 300;
b = 30;
cbar = 11.32;
xcgr = 0.35;
hx = 160;
%%hx = 0;
rad2deg = 180/pi;

% Assign state & control variables
vt = x(1);
alpha = x(2)*rad2deg;
beta = x(3)*rad2deg;
phi = x(4);
theta = x(5);
psi = x(6);
p = x(7);
q = x(8);
r = x(9);
alt = x(12);
pow = x(13);

% air data computer and engine model
[amach,qbar] = adcf16(vt,alt);
cpow = tgear_f16(thtl);
xd(13,1) = pdot_f16(pow,cpow);
t = thrustf16(pow,alt,amach);

% look-up tables and component build-up
cxt = cxf16(alpha,el);
cyt = cyf16(beta,ail,rdr);
czt = czf16(alpha,beta,el);
dail = ail/20;
drdr = rdr/30;
clt = clf16(alpha,beta) + dldaf16(alpha,beta)*dail ...
    + dldrf16(alpha,beta)*drdr;
cmt = cmf16(alpha,el);
cnt = cnf16(alpha,beta) + dndaf16(alpha,beta)*dail ...
    + dndrf16(alpha,beta)*drdr;

% add damping derivatives
tvt = 0.5/vt;
b2v = b*tvt;
cq = cbar*q*tvt;
d = dampf16(alpha);
cxt = cxt + cq*d(1);
cyt = cyt + b2v*(d(2)*r + d(3)*p);
czt = czt + cq*d(4);
clt = clt + b2v*(d(5)*r + d(6)*p);
cmt = cmt + cq*d(7) + czt*(xcgr - xcg);
cnt = cnt + b2v*(d(8)*r + d(9)*p) - cyt*(xcgr - xcg)*cbar/b;

% get ready for state equations
cbta = cos(x(3));
u = vt*cos(x(2))*cbta;
v = vt*sin(x(3));
w = vt*sin(x(2))*cbta;

sth = sin(theta);
cth = cos(theta);
sph = sin(phi);
cph = cos(phi);
spsi = sin(psi);
cpsi = cos(psi);
qs = qbar*s;
qsb = qs*b;
rmqs = qs/mass;
gcth = gd*cth;
qsph = q*sph;
ay = rmqs*cyt;
az = rmqs*czt;

% force equations
udot = r*v - q*w - gd*sth + (qs*cxt + t)/mass;
vdot = p*w - r*u + gcth*sph + ay;
wdot = q*u - p*v + gcth*cph + az;
dum = u*u + w*w;
xd(1,1) = (u*udot + v*vdot + w*wdot)/vt;
xd(2,1) = (u*wdot - w*udot)/dum;
xd(3,1) = (vt*vdot - v*xd(1))*cbta/dum;

% kinematics
xd(4,1) = p + (sth/cth)*(qsph + r*cph);
xd(5,1) = q*cph - r*sph;
xd(6,1) = (qsph + r*cph)/cth;

% moments
roll = qsb*clt;
pitch = qs*cbar*cmt;
yaw = qsb*cnt;
pq = p*q;
qr = q*r;
qhx = q*hx;
xd(7,1) = (xpq*pq - xqr*qr + azz*roll + axz*(yaw + qhx))/gam;
xd(8,1) = (ypr*p*r - axz*(p*p - r*r) + pitch - r*hx)/ayy;
xd(9,1) = (zpq*pq - xpq*qr + axz*roll + axx*(yaw + qhx))/gam;

% navigation
t1 = sph*cpsi;
t2 = cph*sth;
t3 = sph*spsi;
smat = [cth*cpsi  t1*sth-cph*spsi  t2*cpsi+t3;
        cth*spsi  t3*sth+cph*cpsi  t2*spsi-t1;
        sth       -sph*cth         -cph*cth];

xd(10:12,1) = smat*[u; v; w];
            % xd(10) = north speed
            % xd(11) = east speed
            % xd(12) = vertical speed

% other outputs
an = -az/gd;
alat = ay/gd;
