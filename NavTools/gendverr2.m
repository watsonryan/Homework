function dverr = gendverr2(deltav,time,errparam,dvseed)
%GENDVERR2		Delta-V error generator.  Note that the error parameters
%                       are defined differently than in GENDVERR.M
%
%	dverr = gendverr2(deltav,time,errparam,dvseed)
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
%                  errparam(2,1) = x accel scale factor error (in ppm)
%                  errparam(2,2) = y accel scale factor error (in ppm)
%                  errparam(2,3) = z accel scale factor error (in ppm)
%                  errparam(3,1) = x accel velocity random walk (max)
%                                  (in meters per second per root-hour)
%                  errparam(3,2) = y accel velocity random walk (max)
%                                  (in meters per second per root-hour)
%                  errparam(3,3) = z accel velocity random walk (max)
%                                  (in meters per second per root-hour)
%
%       dvseed = Optional seed for Gaussian random number generator.
%                Default setting is: saseed = sum(100*clock).
%
%   OUTPUTS
%       dverr = profile of delta-V errors (in body frame: Nose-Rt.Wing-Down)
%

%	M. & S. Braasch 7-98;  Revised 09-2010
%	Copyright (c) 1998-2010 by GPSoft
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

mpsprh2mpsprs = 1/sqrt(3600);  % conversion factor going from
%                          % meters-per-second-per-root-hour to
%                          % meters-per-second-per-root-second

xbias = tdvec*errparam(1,1);
ybias = tdvec*errparam(1,2);
zbias = tdvec*errparam(1,3);
xsferr = deltav(:,1)*errparam(2,1)*1e-6;
ysferr = deltav(:,2)*errparam(2,2)*1e-6;
zsferr = deltav(:,3)*errparam(2,3)*1e-6;
xnoise = ( sqrt(tdvec).*2.*(rand(K,1)-0.5) )*errparam(3,1)*mpsprh2mpsprs;
ynoise = ( sqrt(tdvec).*2.*(rand(K,1)-0.5) )*errparam(3,2)*mpsprh2mpsprs;
znoise = ( sqrt(tdvec).*2.*(rand(K,1)-0.5) )*errparam(3,3)*mpsprh2mpsprs;

dverr = [xbias+xsferr+xnoise ybias+ysferr+ynoise zbias+zsferr+znoise];
