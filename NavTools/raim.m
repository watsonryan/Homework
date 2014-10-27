function parvec = raim(prvec,svxyzmat,initpos,tol)
%RAIM	Receiver Autonomous Integrity Monitoring.  RAIM computes the parity vector 
%       associated with the overdetermined user position solution.
%       
%	parvec = RAIM(prvec,svxyzmat,initpos,tol)
%
%   INPUTS
%	prvec = vector of 'measured' pseudoranges for satellites
%               specified in svxyzmat
%	svxyzmat(i,1:3) = position of satellite i in user defined 
%                         cartesian coordinates.
%	initpos = optional argument.  Initial 'estimate' of user state:
%                 three-dimensional position and clock offset
%                 (in user defined coordinates).  Used to speed up 
%                 iterative solution.  Initial clock offset is
%                 optional with default value = 0.
%	tol = optional argument.  Tolerance value used to determine 
%             convergence of iterative solution.  Default value = 1e-3
%
%   OUTPUTS
%	parvec = parity vector with n-4 elements where n is the total number
%                of pseudorange measurements

%	References: 
%                   Understanding GPS: Principles and Applications,
%	            Elliott D. Kaplan, Editor, Artech House Publishers,
%	            Boston, 1996.
%
%	M. & S. Braasch 12-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.
%

if nargin<4,tol=1e-3;end
if nargin<3,initpos=[0 0 0 0];end
if nargin<2,error('insufficient number of input arguments'),end
[m,n]=size(initpos);
if m>n, estusr=initpos';else,estusr=initpos;end
if max(size(estusr))<3,
   error('must define at least 3 dimensions in INITPOS')
end
if max(size(estusr))<4,estusr=[estusr 0];end
numvis=max(size(svxyzmat));
beta=[1e9 1e9 1e9 1e9];
maxiter=10;
iter=0;
while ((iter<maxiter)&(norm(beta)>tol)),
    for N = 1:numvis,
	pr0 = norm(svxyzmat(N,:)-estusr(1:3));
	y(N,1) = prvec(N) - pr0 - estusr(4);
    end,
    H = hmat(svxyzmat,estusr(1:3));
    beta = H\y;
    estusr=estusr+beta';
    iter=iter+1;
end
      [Q,R] = qr(H);
      Qt = Q';
%%      Q_beta = Qt(1:3,:);
%%      U = R(1:3,:);
      [r,c] = size(Qt);
      Q_parity = Qt(5:r,:);
      parvec=Q_parity*y;

