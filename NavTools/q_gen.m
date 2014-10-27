function [Q,PHI] = q_gen(deltat,F,G,W,N)
%Q_GEN	Compute the Q and PHI matrices for an N-state Kalman filter
%
%  [Q,PHI] = q_gen(deltat,F,G,W,N)
%
%   INPUTS
%	deltat = filter time step (update interval) in seconds
%   F = continuous-time system dynamics matrix
%   G = input coefficient matrix
%   W = input power spectral density matrix
%   N = number of states in the filter
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

if nargin<5,error('insufficient number of input arguments'),end

TwoN = 2*N;
A = zeros(TwoN,TwoN);
A(1:N,1:N) = -1*F;
A(1:N,N+1:TwoN) = G*W*G';
A(N+1:TwoN,N+1:TwoN) = F';
A = A*deltat;

B = expm(A);

PHI_trans = B(N+1:TwoN,N+1:TwoN);
PHI = PHI_trans';

Q = PHI*B(1:15,16:30);
