function rt = rtau(dp)
%RTAU   Reciprocal time constant
%
% Used by the function PDOT_F16
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

if dp <= 25,
    rt = 1;        % reciprocal time constant
elseif dp >= 50,
    rt = 0.1;
else
    rt = 1.9 - 0.036*dp;
end
