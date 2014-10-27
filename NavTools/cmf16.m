function pm_coeff = cmf16(alpha,el)
%CMF16  F-16 pitching moment coefficient
%
%  INPUTS
%     alpha = angle-of-attack in degrees
%     el = elevator deflection angle in degrees
%
%  OUTPUT
%     pm_coeff = pitching moment coefficient
%
%  In the table, ALPHA ranges from -10 to 45 degrees
%  in 5 degree increments.  EL ranges from -24 degrees
%  to +24 degrees in 12 degree increments
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

a = [0.205  0.168  0.186  0.196  0.213  0.251  0.245 ...
     0.238  0.252  0.231  0.198  0.192;...
     0.081  0.077  0.107  0.110  0.110  0.141  0.127 ...
     0.119  0.133  0.108  0.081  0.093;...
     -0.046  -0.02  -0.009  -0.005  -0.006  0.01  0.006 ...
     -0.001  0.014  0.0  -0.013  0.032;...
     -0.174  -0.145  -0.121  -0.127  -0.129  -0.102  -0.097 ...
     -0.113  -0.087  -0.084  -0.069  -0.006;...
     -0.259  -0.202  -0.184  -0.193  -0.199  -0.150  -0.160 ...
     -0.167  -0.104  -0.076  -0.041  -0.005];

  s = 0.2*alpha;
  k = fix(s);
  if k <= -2, 
      k = -1;
  end
  if k >= 9, 
      k = 8;
  end
  da = s - k;
  L = k + sign(da);
  
  s = el/12;
  m = fix(s);
  if m <= -2, m=-1; end
  if m >= 2; m=1; end
  de = s - m;
  n = m + sign(de);
  
  k = k + 3;
  L = L + 3;
  m = m + 3;
  n = n + 3;
  t = a(m,k);
  u = a(n,k);
 
  v = t + abs(da)*(a(m,L)-t);
  w = u + abs(da)*(a(n,L)-u);
  pm_coeff = v + (w-v)*abs(de);
  