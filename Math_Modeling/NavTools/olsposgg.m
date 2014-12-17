function estusr = olsposgg(prvec,svxyzmat,svid,initpos,tol)
%OLSPOSGG	Compute position from satellite positions and pseudoranges
%       	via ordinary least squares.  Routine solves for receiver clock
%               offset and system time offset (e.g., between GPS
%             and Galileo)
%
%	estusr = OLSPOSGG(prvec,svxyzmat,svid,initpos,tol)
%
%   INPUTS
%	prvec = vector of 'measured' pseudoranges for satellites
%               specified in svxyzmat
%	svxyzmat(i,1:3) = position of satellite i in user defined 
%                         cartesian coordinates.
%       svid = vector of satellite identification numbers corresponding to svxyzmat
%	initpos = optional argument.  Initial 'estimate' of user state:
%                 three-dimensional position and clock offset
%                 (in user defined coordinates).  Used to speed up 
%                 iterative solution.  Initial clock offset is
%                 optional with default value = 0.
%	tol = optional argument.  Tolerance value used to determine 
%             convergence of iterative solution.  Default value = 1e-3
%
%   OUTPUTS
%	estusr(1:3) = estimated user x, y, and z coordinates
%	estusr(4) = estimated user clock offset
%       estusr(5) = estimate GPS/Glonass or GPS/Galileo time offset
%          Note: all five elements of estusr are in the same units
%                as those used in prvec

%	M. & S. Braasch 12-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.
%

if nargin<5,tol=1e-3;end
if nargin<4,initpos=[0 0 0 0];end
if nargin<3,error('insufficient number of input arguments'),end
[m,n]=size(initpos);
if m>n, estusr=initpos';else,estusr=initpos;end
if max(size(estusr))<3,
   error('must define at least 3 dimensions in INITPOS')
end
if max(size(estusr))<4,estusr=[estusr 0];end
if max(size(estusr))<5,estusr=[estusr 0];end
numvis=max(size(svxyzmat));
if numvis<5, error('Cannot compute hybrid GPS/Glonass solution with less than 5 satellites'),end
beta=[1e9 1e9 1e9 1e9 1e9];
maxiter=10;
iter=0;
while ((iter<maxiter)&(norm(beta)>tol)),
    for N = 1:numvis,
	pr0 = norm(svxyzmat(N,:)-estusr(1:3));
	y(N,1) = prvec(N) - pr0 - estusr(4);
        if svid(N) > 49, y(N,1) = y(N,1) - estusr(5); end
    end,
    H = hmat(svxyzmat,estusr(1:3));
    hz = zeros(numvis,1);  H = [H hz];
    for i = 1:numvis,
        if svid(i) > 49, H(i,5) = 1; end
    end
    beta = H\y;
    estusr=estusr+beta';
    iter=iter+1;
end

