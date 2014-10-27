function [q,phi] = q_gen_5_gps(deltat,sp,sf,sg)
%Q_GEN_5_GPS	Compute the Q and PHI matrices for a 5-state
%               GPS stand alone extended Kalman filter
%
%  [q,phi] = q_gen_5_gps(deltat,sp,sf,sg)
%
%   INPUTS
%	deltat = filter time step (update interval) in seconds
%   sp = spectral amplitude of white-noise input to position random-walk process
%   sf = spectral amplitude of clock frequency noise
%   sg = spectral amplitude of clock frequency rate noise
%
%   OUTPUTS
%   q = system noise covariance matrix
%   phi = state transition matrix

%	Reference:
%   Levy, L., "Integration of GPS with Inertial Navigation Systems,"
%             short course notes, Navtech Seminars, Springfield, VA, 2003.
%
%	Copyright (c) 2005 by GPSoft LLC


F=zeros(5,5); F(4,5)=1;

rsp = sqrt(sp);  % square root of the spectral amplitude for the position state noise
rsf = sqrt(sf);  % square root of the spectral amplitude for the clock bias state noise
rsg = sqrt(sg);  % square root of the spectral amplitude for the clock rate state noise

G=zeros(5,5); G(1,1)=rsp; G(2,2)=rsp; G(3,3)=rsp; G(4,4)=rsf; G(5,5)=rsg;
W = zeros(5,5); W(1,1)=1; W(2,2)=1; W(3,3)=1; W(4,4)=1; W(5,5)=1;

A = zeros(10,10);
A(1:5,1:5) = -1*F;
A(1:5,6:10) = G*W*G';
A(6:10,6:10) = F';
A = A*deltat;

B = expm(A);

PHI_trans = B(6:10,6:10);
phi = PHI_trans';

q = phi*B(1:5,6:10);
