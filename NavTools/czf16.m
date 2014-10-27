function force = czf16(alpha,beta,el)
%CZF16  F-16 z-axis force coefficient
%
% force = czf16(alpha,beta,el)
%
%  INPUTS
%     alpha = angle-of-attack in degrees
%     beta = sideslip angle in degrees
%     el = elevator deflection angle in degrees
%
%  OUTPUT
%     force = z-axis force coefficient
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%

a = [.770 .241 -.100 -.416 -.731 -1.053 -1.366 -1.646 ...
        -1.917 -2.120 -2.248 -2.229];

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
  
  s = a(k) + abs(da)*(a(L)-a(k));
  force = s*(1-(beta*pi/180)^2) - 0.19*(el/25);
  