function v_ex = extrapol(v1,v2,td12,tdex)
%EXTRAPOL      Perform extrapolation.  The values v1 and v2
%              are spaced td12 seconds apart in time.  Assume
%              v1 corresponds to a time of t1 and that v2
%              corresponds to a time t2 (t2 > t1).  The
%              function extrapolates these values to a point
%              in time which is tdex seconds after t2.
%
%  v_ex = extrapol(v1,v2,td12,tdex);
%

%	M. & S. Braasch 4-98
%	Copyright (c) 1997-98 by GPSoft
%	All Rights Reserved.
%

if nargin<4,error('insufficient number of input arguments'),end

v_ex = v2 + ( (v2 - v1)/td12 )*tdex;
