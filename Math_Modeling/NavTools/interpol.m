function v_in = interpol(v1,v2,td12,tdin)
%INTERPOL      Perform interpolation.  The values v1 and v2
%              are spaced td12 seconds apart in time.  Assume
%              v1 corresponds to a time of t1 and that v2
%              corresponds to a time t2 (t2 > t1).  The
%              function interpolates these values to a point
%              in time which is tdin seconds after t1.
%
%  v_in = interpol(v1,v2,td12,tdin);
%

%	M. & S. Braasch 6-98
%	Copyright (c) 1998 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<4,error('insufficient number of input arguments'),end

v_in = v1 + ( (v2 - v1)/td12 )*tdin;
