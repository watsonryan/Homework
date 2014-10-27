function cost = costf16(s)
%COSTF16    Scalar cost function for the F-16 6DOF model
%
%  cost = costf16(s)
%
%  INPUT
%    s = parameter vector
%        s(1) = throttle setting (ranges from 0 to 1)
%        s(2) = elevator deflection in degrees
%        s(3) = angle of attack in radians
%        s(4) = aileron deflection in degrees
%        s(5) = rudder deflection in degrees
%        s(6) = sideslip angle in radians
%
%  OUTPUT
%    cost = scalar cost function
%        cost is a weighted sum of squares of certain elements in the
%        derivative of the state vector.  These elements are: airspeed,
%        angle of attack, sideslip angle, roll-rate, pitch-rate and
%        yaw-rate.
%
%  GLOBAL VARIABLES
%    See the program for details.  The global variables include the
%    aircraft state vector and its derivative, climb angle, roll-rate,
%    pitch-rate, turn-rate, turn coordination flag.
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%
%    
global xcg x xd gamma rollrate pitchrate turnrate coord stab
global an alat qbar amach q alpha

thtl = s(1);
el = s(2);
x(2) = s(3);
ail = s(4);
rdr = s(5);
x(3) = s(6);
x(13) = tgear_f16(thtl);
x = constr_f16(x,gamma,rollrate,pitchrate,turnrate,coord,stab);
time = 0;
[xd,an,alat,qbar,amach,q,alpha] = ...
         f16_6dof(time,x,thtl,el,ail,rdr,xcg);
cost = xd(1)^2 + 100*(xd(2)^2 + xd(3)^2) + ...
    10*(xd(7)^2 + xd(8)^2 + xd(9)^2);

