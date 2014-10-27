function qua_vec = dcm2qua(DCMbn)
%DCM2QUA       Direction cosine matrix to quaternion conversion.
%       
%	qua_vec = dcm2qua(DCMbn)
%
%   INPUTS
%       DCMbn = 3x3 direction cosine matrix providing the
%             transformation from the body frame
%             to the navigation frame
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
%   NOTE
%       The algorithm assumes small angular displacements.
%

%   REFERENCE:  Titterton, D. and J. Weston, STRAPDOWN
%               INERTIAL NAVIGATION TECHNOLOGY, Peter
%               Peregrinus Ltd. on behalf of the Institution
%               of Electrical Engineers, London, 1997.
%
%	M. & S. Braasch 12-97
%	Copyright (c) 1997 by GPSoft
%	All Rights Reserved.
%

  if nargin<1,error('insufficient number of input arguments'),end
  
  a = 0.5*sqrt(1+DCMbn(1,1)+DCMbn(2,2)+DCMbn(3,3));
  tmp = inv(4*a);
  b = tmp*(DCMbn(3,2)-DCMbn(2,3));
  c = tmp*(DCMbn(1,3)-DCMbn(3,1));
  d = tmp*(DCMbn(2,1)-DCMbn(1,2));
  
  qua_vec = [a b c d];
  