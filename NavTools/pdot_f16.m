function pd = pdot_f16(p3,p1)
%PDOT_F16     Compute rate-of-change of power
%
%  pd = pdot_f16(p3,p1)
%
%  INPUTS
%     p3 = actual power
%     p1 = commanded power
%
%  OUTPUT
%     pd = rate-of-change of power
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

if p1 >= 50,
    if p3 >= 50,
        t = 5;
        p2 = p1;
    else
        p2 = 60;
        t = rtau(p2-p3);
    end
else
    if p3 >= 50,
        t = 5;
        p2 = 40;
    else
        p2 = p1;
        t = rtau(p2-p3);
    end
end

pd = t*(p2-p3);
