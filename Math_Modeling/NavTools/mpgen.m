function delta=mpgen(numsat,duration,model,mpseed)
%MPGEN	Non-"Realtime" generation of zero-angle equivalent
%       pseudorange multipath error
%
%       delta = MPGEN(numsat,duration,model,mpseed)
%
%   INPUTS
%       numsat = number of multipath error traces (i.e., number of 
%                satellites for which multipath is to be generated)
%	duration = duration of error traces in seconds (must be an integer)
%       model = Optional model choice for generation of 
%               multipath.  Default is: model = 1
%               1 = standard deviation is approximately 1.6 meters
%		    time constant is approximately 2 minutes
%               2 = standard deviation is approximately 0.40 meters
%		    time constant is approximately 15 seconds
%	mpseed = Optional seed for Gaussian random number generator.
%                Default setting is: mpseed = sum(100*clock).
%
%   OUTPUTS
%       delta = zero-angle equivalent multipath error in meters
%
%   NOTE:  All error traces are generated at a rate of 1 Hz.

%
%	M. & S. Braasch 11-96;  Revised 10-99
%	Copyright (c) 1996-99 by GPSoft LLC
%	All Rights Reserved.
%
if nargin<4,randn('seed',sum(100*clock)),else,randn('seed',mpseed);end
if nargin<3,model=1;end
if nargin<2,error('insufficient number of input arguments'),end

% Long time constant Autoregressive Model
   if model==1,
      %[b,a] = butter(1,0.007);
      b = [0.01087642013487   0.01087642013487];
      a = [1.00000000000000  -0.97824715973025];
      sigmae = 15;
      for k = 1:numsat,
          x = sigmae*randn(duration,1);
          y = filter(b,a,x);
          delta(:,k) = y;
      end
   end,
%
% Short time constant Autoregressive Model
   if model==2,
      %[b,a] = butter(1,0.025);
      b = [0.03780475417090   0.03780475417090];
      a = [1.00000000000000  -0.92439049165821];
      sigmae = 2;
      for k = 1:numsat,
          x = sigmae*randn(duration,1);
          y = filter(b,a,x);
          delta(:,k) = y;
      end
   end,

