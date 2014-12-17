function [xnew] = rk4_f16(time,dt,xx,thtl,el,ail,rdr,xcg)
%RK4_F16    Fourth-order Runge-Kutta integration algorithm to solve the
%           F-16 6DOF model and update the state vector
%
%  xnew = rk4_f16(time,dt,xx,thtl,el,ail,rdr,xcg)
%
%  INPUTS
%     time = time in seconds
%     dt = integration time step in seconds
%     xx = aircraft state vector
%     thtl = throttle setting (ranges from 0 to 1)
%     el = elevator deflection in degrees
%     ail = aileron deflection in degrees
%     rdr = rudder deflection in degrees
%     xcg = center of gravity position relative to the leading edge of the
%           wing in terms of the mean aerodynamic chord (mac); nominal is
%           0.35mac; a forward cg is 0.3mac and an aft cg is 0.38mac
%
%  OUTPUT
%     xnew = updated aircraft state vector
%
%  REFERENCE
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

[xd] = f16_6dof(time,xx,thtl,el,ail,rdr,xcg);

xa = xd*dt;
x = xx + 0.5*xa;
t = time + 0.5*dt;
[xd] = f16_6dof(t,x,thtl,el,ail,rdr,xcg);

q = xd*dt;
x = xx + 0.5*q;
xa = xa + 2*q;
[xd] = f16_6dof(t,x,thtl,el,ail,rdr,xcg);

q = xd*dt;
x = xx + q;
xa = xa + 2*q;
time = time + dt;
[xd] = f16_6dof(time,x,thtl,el,ail,rdr,xcg);

xnew = xx + (xa + xd*dt)/6;
