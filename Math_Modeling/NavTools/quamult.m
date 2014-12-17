function quanew = quamult(qua1,qua2)
%QUAMULT       Quaternion multiplication.
%       
%	 quanew = quamult(qua1,qua2)
%
%   INPUT
%       qua1, qua2 = 4 element quaternion vectors
%                  = [a b c d]
%       where: a = cos(MU/2)
%              b = (MUx/MU)*sin(MU/2)
%              c = (MUy/MU)*sin(MU/2)
%              d = (MUz/MU)*sin(MU/2)
%       where: MUx, MUy, MUz are the components of the angle vector
%              MU is the magnitude of the angle vector
%
%   OUTPUT
%       quanew = (qua1)*(qua2) using quaternion multiplication
%

%   REFERENCE:  Titterton, D. and J. Weston, STRAPDOWN
%               INERTIAL NAVIGATION TECHNOLOGY, Peter
%               Peregrinus Ltd. on behalf of the Institution
%               of Electrical Engineers, London, 1997.
%
%	M. & S. Braasch 1-98
%	Copyright (c) 1998 by GPSoft
%	All Rights Reserved.
%

  if nargin<2,error('insufficient number of input arguments'),end
  [m,n]=size(qua2); if m>n, quatmp=qua2; else, quatmp=qua2'; end
  
  a = qua1(1); b = qua1(2); c = qua1(3); d = qua1(4);
  
  A = [ a -b -c -d;
        b  a -d  c;
        c  d  a -b;
        d -c  b  a];
    
  quanew = (A*quatmp)';
  