function [Q,PHI] = q_gen_8_gps(deltat,Sp,Sf,Sg)
%Q_GEN_8_GPS	Compute the Q and PHI matrices for a 8-state
%               GPS stand alone extended Kalman filter
%
%  [Q,PHI] = q_gen_8_gps(deltat,sp,sf,sg)
%
%   INPUTS
%	deltat = filter time step (update interval) in seconds
%   Sp = spectral amplitude of white-noise input to position random-walk process
%   Sf = spectral amplitude of clock frequency noise
%   Sg = spectral amplitude of clock frequency rate noise
%
%   OUTPUTS
%   q = system noise covariance matrix
%   phi = state transition matrix

%	Reference: Brown, R. G. and P. Y. C. Hwang, "Introduction to Random
%	Signals and Applied Kalman Filtering," 3rd edition, John Wiley & Sons,
%	1997.
%
%	Copyright (c) 2005 by GPSoft LLC


F = zeros(8,8); F(1,2)=1; F(3,4)=1; F(5,6)=1; F(7,8)=1;
PHI = expm(F*deltat);

Q = zeros(8,8); 
Qsub1 = [Sp*(deltat^3)/3 Sp*(deltat^2)/2;
        Sp*(deltat^2)/2 Sp*deltat];
Qsub2 = [Sf*deltat + Sg*(deltat^3)/3  Sg*(deltat^2)/2;
         Sg*(deltat^2)/2                Sg*deltat];
Q(1:2,1:2) = Qsub1;
Q(3:4,3:4) = Qsub1;
Q(5:6,5:6) = Qsub1;
Q(7:8,7:8) = Qsub2;
