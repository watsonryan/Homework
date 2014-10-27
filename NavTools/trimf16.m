function [init_cost,final_cost] = trimf16(thtl,el,ail,rdr)
%TRIMF16    Function to solve for trim conditions in the F-16 6DOF model
%
%  [init_cost,final_cost] = trimf16(thtl,el,ail,rdr)
%
%  INPUTS
%     thtl = throttle setting (ranges from 0 to 1)
%     el = elevator deflection in degrees
%     ail = aileron deflection in degrees
%     rdr = rudder deflection in degrees
%
%  OUTPUT
%     init_cost = initial value of cost function prior to optimization
%     final_cost = final value of cost function after optimization
%
%  GLOBAL VARIABLES
%     See program for details.  The global variables include the aircraft
%     state vector, its derivative, roll-rate, pitch-rate, turn-rate, turn
%     coordination flag, and final positions of the throttle, elevator,
%     aileron and rudder.
%
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
global thtl_final el_final ail_final rdr_final

s0(1) = thtl;
s0(2) = el;
s0(3) = x(2);
s0(4) = ail;
s0(5) = rdr;
s0(6) = x(3);
init_cost = costf16(s0);
options = optimset('MaxFunEvals',10000,'MaxIter',10000,'TolFun',1e-12);

for i = 1:5,
    [s,final_cost] = fminsearch('costf16',s0,options);
    s0 = s;
end

thtl_final = s(1);
el_final = s(2);
ail_final = s(4);
rdr_final = s(5);
