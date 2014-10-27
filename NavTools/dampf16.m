function d = dampf16(alpha)
%DAMPF16    Function to compute various damping coefficients
% 
%  d = dampf16(alpha)
%
%  INPUT
%    alpha = angle of attack in degrees
%
%  OUTPUT
%    d = vector of nine damping coefficients
%        d(1) = CXq
%        d(2) = CYr
%        d(3) = CYp
%        d(4) = CZq
%        d(5) = Clr
%        d(6) = Clp
%        d(7) = Cmq
%        d(8) = Cnr
%        d(9) = Cnp
%
%  In the look-up table, ALPHA ranges from
%  -10 to 45 degrees in increments of 5 degrees.
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

a = [-.267  -.110  .308  1.34  2.08  2.91  2.76 ...
      2.05  1.50  1.49 1.83  1.21; ...
      .882  .852  .876  .958  .962  .974  .819 ...
      .483  .590  1.21  -.493  1.04;...
      -.108  -.108  -.188  .110  .258  .226 .344 ...
      .362  .611  .529  .298  -2.27;...
      -8.8  -25.8  -28.9  -31.4  -31.2  -30.7  -27.7 ...
      -28.2  -29.0  -29.8  -38.3  -35.3;...
      -.126  -.026  .063  .113  .208  .230  .319 ...
      .437  .680  .100  .447  -.330;...
      -.360  -.359  -.443  -.420  -.383  -.375  -.329 ...
      -.294  -.230  -.210  -.120  -.100;...
      -7.21  -.540  -5.23  -5.26  -6.11  -6.64  -5.69...
      -6  -6.2  -6.4  -6.6  -6;...
      -.380  -.363  -.378  -.386  -.370  -.453  -.550 ...
      -.582  -.595  -.637  -1.02  -.840;...
      .061  .052  .052  -.012  -.013  -.024  .050...
      .150  .130  .158  .240  .150];

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
  
  k = k + 3;
  L = L + 3;
  d = a(:,k) + abs(da)*( a(:,L) - a(:,k) );
  