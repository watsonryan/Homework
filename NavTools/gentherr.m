function dtherr = gentherr(deltath,time,errparam,dthseed)
%GENTHERR		Delta-theta error generator. 
%       
%	dtherr = gentherr(deltath,time,errparam,dthseed)
%
%   INPUTS
%       deltath = profile of ideal rate-integrating gyro 
%                 outputs over time
%         deltath(i,1:3) = for the i-th flight path segment,
%                          deltath(i,1) is the x gyro data
%                          (nose positive); deltath(i,2) is the y
%                          gyro data (right wing positive);
%                          deltath(i,3) is the z gyro data
%                          (down is positive)
%
%       time = sequential time vector for the simulated flight path (seconds)
%
%       errparam = error parameter matrix
%                  errparam(1,1) = x gyro bias (in degrees per hour)
%                  errparam(1,2) = y gyro bias (in degrees per hour)
%                  errparam(1,3) = z gyro bias (in degrees per hour)
%                  errparam(2,1) = x gyro scale factor error (in percent)
%                  errparam(2,2) = y gyro scale factor error (in percent)
%                  errparam(2,3) = z gyro scale factor error (in percent)
%                  errparam(3,1) = x gyro white noise standard deviation
%                                  (in degrees per square-root of hour)
%                  errparam(3,2) = y gyro white noise standard deviation
%                                  (in degrees per square-root of hour)
%                  errparam(3,3) = z gyro white noise standard deviation
%                                  (in degrees per square-root of hour)
%
%       dthseed = Optional seed for Gaussian random number generator.
%                Default setting is: saseed = sum(100*clock).
%
%   OUTPUTS
%       dtherr = profile of delta-theta errors 
%                (in body frame: Nose-Rt.Wing-Down)
%

%	M. & S. Braasch 7-98
%	Copyright (c) 1998 by GPSoft
%	All Rights Reserved.
%

if nargin<4,randn('seed',sum(100*clock)),else,randn('seed',dthseed);end
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

dph2rps = (pi/180)/3600;   % conversion constant from deg/hr to rad/sec

dprh2rprs = (pi/180)/sqrt(3600);  % conversion factor going from
%                                 % degrees-per-root-hour to
%                                 % radians-per-root-second

xbias = tdvec*errparam(1,1)*dph2rps;
ybias = tdvec*errparam(1,2)*dph2rps;
zbias = tdvec*errparam(1,3)*dph2rps;
xsferr = deltath(:,1)*errparam(2,1)*0.01;
ysferr = deltath(:,2)*errparam(2,2)*0.01;
zsferr = deltath(:,3)*errparam(2,3)*0.01;
xnoise = ( sqrt(tdvec).*randn(K,1) )*errparam(3,1)*dprh2rprs;
ynoise = ( sqrt(tdvec).*randn(K,1) )*errparam(3,2)*dprh2rprs;
znoise = ( sqrt(tdvec).*randn(K,1) )*errparam(3,3)*dprh2rprs;

dtherr = [xbias+xsferr+xnoise ybias+ysferr+ynoise zbias+zsferr+znoise];
