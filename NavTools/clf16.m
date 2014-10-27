function rm_coeff = clf16(alpha,beta)
%CLF16  F-16 rolling moment coefficient
%
%  rm_coeff = clf16(alpha,beta)
%
%  INPUTS
%     alpha = angle-of-attack in degrees
%     beta = sideslip angle in degrees
%
%  OUTPUT
%     rm_coeff = rolling moment coefficient
%
% In the look-up table, ALPHA ranges from -10 to 45 degrees
% in increments of 5 degrees.  ABS(BETA) ranges from 0 to 30 degrees
% in increments of 5 degrees
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

a = [   0       0        0      0       0       0      0   ...
        0       0        0      0       0;...
     -0.001  -0.004  -0.008  -0.012  -0.016  -0.019  -0.02 ...
      -0.02  -0.015  -0.008  -0.013  -0.015;...
      -0.003  -0.009  -0.017  -0.024  -0.03  -0.034  -0.04 ...
      -0.037  -0.016  -0.002  -0.01  -0.019;...
      -0.001  -0.01  -0.02  -0.03  -0.039  -0.044  -0.05 ...
      -0.049  -0.023  -0.006  -0.014  -0.027;...
      0.00  -0.01  -0.022  -0.034  -0.047  -0.046  -0.059 ...
      -0.061  -0.033  -0.036  -0.035  -0.035;...
       0.007  -0.01  -0.023  -0.034  -0.049  -0.046  -0.068 ...
       -0.071  -0.06  -0.058  -0.062  -0.059;...
       -0.009  -0.011  -0.023  -0.037  -0.05  -0.047  -0.074 ...
       -0.079  -0.091  -0.076  -0.077  -0.076];

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
  
  s = 0.2*abs(beta);
  m = fix(s);
  if m == 0, m=1; end
  if m >= 6; m=5; end
  db = s - m;
  n = m + sign(db);
  
  k = k + 3;
  L = L + 3;
  m = m + 1;
  n = n + 1;
  t = a(m,k);
  u = a(n,k);
 
  v = t + abs(da)*(a(m,L)-t);
  w = u + abs(da)*(a(n,L)-u);
  dum = v + (w-v)*abs(db);
  %if beta == 0,
  %    rm_coeff = dum + 1;
  %else
  %    rm_coeff = dum + sign(beta);
  %end
  rm_coeff = dum*sign(beta);
  