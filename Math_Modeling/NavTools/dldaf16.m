function rma = dldaf16(alpha,beta)
%DLDAF16    Rolling moment due to ailerons
%
%  rma = dldaf16(alpha,beta)
%
%  INPUTS
%     alpha = angle-of-attack in degrees
%     beta = sideslip angle in degrees
%
%  OUTPUT
%     rma = rolling moment due to ailerons
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

a = [ -0.041  -0.052  -0.053  -0.056  -0.05  -0.056  -0.082 ...
      -0.059  -0.042  -0.038  -0.027  -0.017;...
      -0.041  -0.053  -0.053  -0.053  -0.05  -0.051  -0.066 ...
      -0.043  -0.038  -0.027  -0.023  -0.016;...
      -0.042  -0.053  -0.052  -0.051  -0.049  -0.049  -0.043 ...
      -0.035  -0.026  -0.016  -0.018  -0.014;...
      -0.04  -0.052  -0.051  -0.052  -0.048  -0.048  -0.042 ...
      -0.037  -0.031  -0.026  -0.017  -0.012;...
      -0.043  -0.049  -0.048  -0.049  -0.043  -0.042  -0.042 ...
      -0.036  -0.025  -0.021  -0.016  -0.011;...
      -0.044  -0.048  -0.048  -0.047  -0.042  -0.041  -0.02 ...
      -0.028  -0.013  -0.014  -0.011  -0.01;...
      -0.043  -0.049  -0.047  -0.045  -0.042  -0.037  -0.003 ...
      -0.013  -0.01  -0.003  -0.007  -0.008];

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
  rma = v + (w-v)*abs(db);
 