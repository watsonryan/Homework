function ym_coeff = cnf16(alpha,beta)
%CNF16  F-16 yawing moment coefficient
%
%  INPUTS
%     alpha = angle-of-attack in degrees
%     beta = sideslip angle in degrees
%
%  OUTPUT
%     ym_coeff = yawing moment coefficient
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

a = [   0      0      0      0      0      0      0   ...
        0      0      0      0      0; ...
      0.018  0.019  0.018  0.019  0.019  0.018  0.013 ...
      0.007  0.004  -0.014  -0.017  -0.033;...
      0.038  0.042  0.042  0.042  0.043  0.039  0.03 ...
      0.017 0.004  -0.035  -0.047  -0.057;...
      0.056  0.057  0.059  0.058  0.058  0.053  0.032 ...
      0.012  0.002  -0.046  -0.071  -0.073;...
      0.064  0.077  0.076  0.074  0.073  0.057  0.029 ...
      0.007  0.012  -0.034  -0.065  -0.041;...
      0.074  0.086  0.093  0.089  0.08  0.062  0.049 ...
      0.022  0.028  -0.012  -0.002  -0.013;...
      0.079  0.09  0.106  0.106  0.096  0.08  0.068 ...
      0.03  0.064  0.015  0.011  -0.001];

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
  %    ym_coeff = dum + 1;
  %else
  %    ym_coeff = dum + sign(beta);
  %end
  ym_coeff = dum*sign(beta);
  