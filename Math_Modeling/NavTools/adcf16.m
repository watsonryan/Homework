function [amach,qbar,ps] = adcf16(vt,alt)
%ADCF16     Standard atmosphere air data computer model for an F-16
%
%   [amach,qbar,ps] = adcf16(vt,alt)
%
%  INPUTS
%     vt = total vehicle velocity in feet/sec
%     alt = altitude in feet
%
%  OUTPUTS
%     amach = mach number
%     qbar = dynamic pressure
%     ps = static pressure
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch Feb-March 2005
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

r0 = 2.377e-3;    % sea-level density
tfac = 1 - alt*0.703e-5;
t = 519*tfac;    % temperature
if alt >= 35000, t = 390; end
rho = r0*(tfac^4.14);    % density
amach = vt/sqrt(1.4*1716.3*t);   % mach number
qbar = 0.5*rho*vt*vt;    % dynamic pressure
ps = 1715*rho*t;   % static pressure
