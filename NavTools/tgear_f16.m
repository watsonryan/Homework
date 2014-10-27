function tg = tgear_f16(thtl)
%TGEAR_F16    Power command versus throttle relationship
%
%  tg = tgear_f16(thtl)
%
%  INPUT
%    thtl = throttle setting (ranges from 0 to 1)
%
%  OUTPUT
%    tg = power command
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

if thtl <= 0.77,
    tg = 64.94*thtl;
else
    tg = 217.38*thtl - 117.38;
end
