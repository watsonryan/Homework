function [Q,PHI] = q_gen_11_gps(deltat,sp,sf,sg,tau_accel)
%Q_GEN_11_GPS	Compute the Q and PHI matrices for an 11-state
%               GPS stand alone Kalman filter
%
%  [Q,PHI] = q_gen_11_gps(deltat,sp,sf,sg)
%
%   INPUTS
%	deltat = filter time step (update interval) in seconds
%   sp = spectral amplitude of position random process
%   sf = spectral amplitude of clock frequency noise
%   sg = spectral amplitude of clock frequency rate noise
%   tau_accel = time constant for first order Gauss-Markov 
%              process describing the acceleration state
%
%   OUTPUTS
%   Q = system noise covariance matrix
%   PHI = state transition matrix

%	Reference:
%   Brown, R. G. and P. Y. C. Hwang, "Introduction to Random Signals and
%   Applied Kalman Filtering," 3rd edition, John Wiley & Sons, New York,
%   1997.
%
%	Copyright (c) 2004-2005 by GPSoft LLC

beta = 1/tau_accel;
F = zeros(11,11); 
Fsub1 = [0 1 0;
         0 0 1;
         0 0 -beta];
Fsub2 = [0 1;
         0 0];
F(1:3,1:3) = Fsub1; F(4:6,4:6)=Fsub1; F(7:9,7:9)=Fsub1; F(10:11,10:11)=Fsub2;

rsp = sqrt(sp);  % square root of the spectral amplitude for the position state noise
rsf = sqrt(sf);  % square root of the spectral amplitude for the clock bias state noise
rsg = sqrt(sg);  % square root of the spectral amplitude for the clock rate state noise

G=zeros(11,11); G(3,3)=rsp; G(6,6)=rsp; G(9,9)=rsp; G(10,10)=rsf; G(11,11)=rsg;
W = zeros(11,11); W(3,3)=1; W(6,6)=1; W(9,9)=1; W(10,10)=1; W(11,11)=1;

A = zeros(22,22);
A(1:11,1:11) = -1*F;
A(1:11,12:22) = G*W*G';
A(12:22,12:22) = F';
A = A*deltat;

B = expm(A);

PHI_trans = B(12:22,12:22);
PHI = PHI_trans';

Q = PHI*B(1:11,12:22);
