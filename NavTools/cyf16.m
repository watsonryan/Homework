function force = cyf16(beta,ail,rdr)
%CYF16  F-16 sideforce coefficient
%
% force = cyf16(beta,ail,rdr)
%
%  INPUTS
%     beta = sideslip angle in degrees
%     ail = aileron deflection angle in degrees
%     rdr = rudder deflection angle in degrees
%
%  OUTPUT
%     force = sideforce coefficient
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

force = -0.02*beta + 0.021*(ail/20) + 0.086*(rdr/30);
