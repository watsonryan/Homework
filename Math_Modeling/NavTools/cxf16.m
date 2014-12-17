function force = cxf16(alpha,el)
%CXF16  F-16 x-axis aerodynamic force coefficient
%
% force = cxf16(alpha,el)
%
%  INPUTS
%     alpha = angle-of-attack in degrees
%     el = elevator deflection angle in degrees
%
%  OUTPUT
%     force = x-axis aerodynamic force coefficient
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

a = [-0.099  -0.081  -0.081  -0.063  -0.025  0.044  0.097 ...
     0.113  0.145  0.167  0.174  0.166;...
     -0.048  -0.038  -0.04  -0.021  0.16  0.083  0.127 ...
     0.137  0.162  0.177  0.179  0.167;...
     -0.022  -0.02  -0.021  -0.004  0.032  0.094  0.128 ...
     0.130  0.154  0.161  0.155  0.138;...
     -0.04  -0.038  -0.039  -0.025  0.006  0.062  0.087 ...
     0.085  0.1  0.11  0.104  0.091;...
     -0.083  -0.073 -0.076  -0.072  -0.046  0.012  0.024 ...
     0.025  0.043  0.053  0.047  0.04];

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
  force = v + (w-v)*abs(de);
  