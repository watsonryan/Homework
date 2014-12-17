function dtherr = gentherr2(deltath,time,errparam,dthseed)
%GENTHERR2		Delta-theta error generator.  Note that the error parameters are defined
%                       differently than in GENTHERR.M
%
%	dtherr = gentherr2(deltath,time,errparam,dthseed)
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
%                  errparam(2,1) = x gyro scale factor error (in ppm)
%                  errparam(2,2) = y gyro scale factor error (in ppm)
%                  errparam(2,3) = z gyro scale factor error (in ppm)
%                  errparam(3,1) = x gyro angle random walk (max)
%                                  (in degrees per root-hour)
%                  errparam(3,2) = y gyro angle random walk (max)
%                                  (in degrees per root-hour)
%                  errparam(3,3) = z gyro angle random walk (max)
%                                  (in degrees per root-hour)
%
%       dthseed = Optional seed for Gaussian random number generator.
%                Default setting is: saseed = sum(100*clock).
%
%   OUTPUTS
%       dtherr = profile of delta-theta errors 
%                (in body frame: Nose-Rt.Wing-Down)
%

%	M. & S. Braasch 7-98;  Revised 09-2010
%	Copyright (c) 1998-2010 by GPSoft
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
xsferr = deltath(:,1)*errparam(2,1)*1e-6;
ysferr = deltath(:,2)*errparam(2,2)*1e-6;
zsferr = deltath(:,3)*errparam(2,3)*1e-6;
xnoise = ( sqrt(tdvec).*2.*(rand(K,1)-0.5) )*errparam(3,1)*dprh2rprs;
ynoise = ( sqrt(tdvec).*2.*(rand(K,1)-0.5) )*errparam(3,2)*dprh2rprs;
znoise = ( sqrt(tdvec).*2.*(rand(K,1)-0.5) )*errparam(3,3)*dprh2rprs;

dtherr = [xbias+xsferr+xnoise ybias+ysferr+ynoise zbias+zsferr+znoise];
