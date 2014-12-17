function dverr = gendverr(deltav,time,errparam,dvseed)
%GENDVERR		Delta-V error generator. 
%       
%	dverr = gendverr(deltav,time,errparam,dvseed)
%
%   INPUTS
%       deltav = profile of ideal accelerometer outputs over time
%          deltav(i,1:3) = for the i-th flight path segment,
%                          deltav(i,1) is the x accelerometer data
%                          (nose positive); deltav(i,2) is the y
%                          accelerometer data (right wing positive);
%                          deltav(i,3) is the z accelerometer data
%                          (down is positive)
%
%       time = sequential time vector for the simulated flight path (seconds)
%
%       errparam = error parameter matrix
%                  errparam(1,1) = x accel bias (in meters per seconds-squared)
%                  errparam(1,2) = y accel bias (in meters per seconds-squared)
%                  errparam(1,3) = z accel bias (in meters per seconds-squared)
%                  errparam(2,1) = x accel scale factor error (in percent)
%                  errparam(2,2) = y accel scale factor error (in percent)
%                  errparam(2,3) = z accel scale factor error (in percent)
%                  errparam(3,1) = x accel white noise standard deviation
%                                  (in meters per seconds-squared)
%                  errparam(3,2) = y accel white noise standard deviation
%                                  (in meters per seconds-squared)
%                  errparam(3,3) = z accel white noise standard deviation
%                                  (in meters per seconds-squared)
%
%       dvseed = Optional seed for Gaussian random number generator.
%                Default setting is: saseed = sum(100*clock).
%
%   OUTPUTS
%       dverr = profile of delta-V errors (in body frame: Nose-Rt.Wing-Down)
%

%	M. & S. Braasch 7-98
%	Copyright (c) 1998 by GPSoft
%	All Rights Reserved.
%

if nargin<4,randn('seed',sum(100*clock)),else,randn('seed',dvseed);end
if nargin<3,error('insufficient number of input arguments'),end

tdint_prof = diff(time);

[m,n] = size(tdint_prof); 
if m<n, 
   tdvec=tdint_prof'; 
   K = n;
else
   tdvec=tdint_prof; 
   K = m;
end

xbias = tdvec*errparam(1,1);
ybias = tdvec*errparam(1,2);
zbias = tdvec*errparam(1,3);
xsferr = deltav(:,1)*errparam(2,1)*0.01;
ysferr = deltav(:,2)*errparam(2,2)*0.01;
zsferr = deltav(:,3)*errparam(2,3)*0.01;
xnoise = ( tdvec.*randn(K,1) )*errparam(3,1);
ynoise = ( tdvec.*randn(K,1) )*errparam(3,2);
znoise = ( tdvec.*randn(K,1) )*errparam(3,3);

dverr = [xbias+xsferr+xnoise ybias+ysferr+ynoise zbias+zsferr+znoise];
