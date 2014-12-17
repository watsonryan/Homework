function ymr = dndrf16(alpha,beta)
%DNDRF16    Yawing moment due to rudder
%
%  ymr = dndrf16(alpha,beta)
%
%  INPUTS
%     alpha = angle-of-attack in degrees
%     beta = sideslip angle in degrees
%
%  OUTPUT
%     ymr = yawing moment due to rudder
%
% In the look-up table, ALPHA ranges from -10 to 45 degrees
% in increments of 5 degrees.  BETA ranges from -30 to 30
% degrees in increments of 10 degrees.
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

a = [-0.018  -0.052  -0.052  -0.052  -0.054  -0.049  -0.059 ...
     -0.051  -0.03  -0.037  -0.026  -0.013;...
     -0.028  -0.051  -0.043  -0.046  -0.045  -0.049  -0.057 ...
     -0.052  -0.03  -0.033  -0.03  -0.008;...
     -0.037  -0.041  -0.038  -0.04  -0.04  -0.038  -0.037 ...
     -0.03  -0.027  -0.024  -0.019  -0.013;...
     -0.048  -0.045  -0.045  -0.045  -0.044  -0.045  -0.047 ...
     -0.048  -0.049  -0.045  -0.033  -0.016;...
     -0.043  -0.044  -0.041  -0.041  -0.04  -0.038  -0.034 ...
     -0.035  -0.035  -0.029  -0.022  -0.009;...
     -0.052  -0.034  -0.036  -0.036  -0.035  -0.028  -0.024 ...
     -0.023  -0.02  -0.016  -0.01  -0.014;...
     -0.062  -0.034  -0.027  -0.028  -0.027  -0.027  -0.023 ...
     -0.023  -0.019  -0.009  -0.025  -0.01];

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
  
  s = 0.1*beta;
  m = fix(s);
  %%if m == -3, m=-2; end
  if m <= -3, m=-2; end
  if m >= 3; m=2; end
  db = s - m;
  n = m + sign(db);
  
  k = k + 3;
  L = L + 3;
  m = m + 4;
  n = n + 4;
  t = a(m,k);
  u = a(n,k);
 
  v = t + abs(da)*(a(m,L)-t);
  w = u + abs(da)*(a(n,L)-u);
  ymr = v + (w-v)*abs(db);
 