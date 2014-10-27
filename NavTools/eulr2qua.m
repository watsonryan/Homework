function qua_vec = eulr2qua(eul_vect)
%EULR2QUA       Euler angle vector to quaternion conversion.
%       
%	qua_vec = eulr2qua(eul_vect)
%
%   INPUTS
%       eul_vect(1) = roll angle in radians 
%
%       eul_vect(2) = pitch angle in radians 
%
%       eul_vect(3) = yaw angle in radians 
%
%   OUTPUTS
%       qua_vec = 4 element quaternion vector
%               = [a b c d]
%       where: a = cos(MU/2)
%              b = (MUx/MU)*sin(MU/2)
%              c = (MUy/MU)*sin(MU/2)
%              d = (MUz/MU)*sin(MU/2)
%       where: MUx, MUy, MUz are the components of the angle vector
%              MU is the magnitude of the angle vector
%

%   October 2009:  Corrected sign error in computation of 'd'
%
%   REFERENCE:  Titterton, D. and J. Weston, STRAPDOWN
%               INERTIAL NAVIGATION TECHNOLOGY, Peter
%               Peregrinus Ltd. on behalf of the Institution
%               of Electrical Engineers, London, 1997.
%
%	M. & S. Braasch December 1997; October 2009
%	Copyright (c) 1997, 2009 by GPSoft
%	All Rights Reserved.
%

  if nargin<1,error('insufficient number of input arguments'),end

  phi = eul_vect(1); theta = eul_vect(2); psi = eul_vect(3);

  cpsi2 = cos(psi/2); spsi2 = sin(psi/2);
  cthe2 = cos(theta/2); sthe2 = sin(theta/2);
  cphi2 = cos(phi/2); sphi2 = sin(phi/2);

  a = cphi2*cthe2*cpsi2 + sphi2*sthe2*spsi2;
  b = sphi2*cthe2*cpsi2 - cphi2*sthe2*spsi2;
  c = cphi2*sthe2*cpsi2 + sphi2*cthe2*spsi2;
  d = cphi2*cthe2*spsi2 - sphi2*sthe2*cpsi2;
  
  qua_vec = [a b c d];
  