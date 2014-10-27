function rmr = dldrf16(alpha,beta)
%DLDRF16    Rolling moment due to rudder
%
%  rmr = dldrf16(alpha,beta)
%
%  INPUTS
%     alpha = angle-of-attack in degrees
%     beta = sideslip angle in degrees
%
%  OUTPUT
%     rmr = rolling moment due to rudder
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
  
  a = [0.005  0.017  0.014  0.01  -0.005  0.009  0.019 ...
       0.005  -0.00  -0.005  -0.011  0.008;...
       0.007  0.016  0.014  0.014  0.013  0.009  0.012 ...
       0.005  0.0  0.004  0.009  0.007;...
       0.013  0.013  0.011  0.012  0.011  0.009  0.008 ...
       0.005  -0.002  0.005  0.003  0.005;...
       0.018  0.015  0.015  0.14  0.014  0.14  0.014 ...
       0.015  0.013  0.011  0.006  0.001;...
       0.015  0.014  0.013  0.013  0.012  0.011  0.011 ...
       0.01  0.008  0.008  0.007  0.003;...
       0.021  0.011  0.01  0.011  0.01  0.009  0.008 ...
       0.01  0.006  0.005  0.0  0.001;...
       0.023  0.01  0.011  0.011  0.011  0.01  0.008 ...
       0.01  0.006  0.014  0.02  0.0];

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
  rmr = v + (w-v)*abs(db);
 