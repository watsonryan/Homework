function deltheta = gendthet(DCMnb_prof)
%GENDTHET      Delta-theta generator.  Local-level components
%              only.  Earth-rate and craft-rate components
%              not calculated (see EARTHROT and CRAFRATE)
%       
%	deltheta = gendthet(DCMnb_prof)
%
%   INPUTS
%       DCMnb_prof = profile of direction cosine elements over time
%                    relating nav-frame to the body-frame
%          DCMnb_prof(i,1:9) = elements of the i-th direction cosine matrix
%                            (DCM) for vehicle attitude; 1 = DCM(1,1),
%                            2 = DCM(1,2), 3 = DCM(1,3),
%                            4 = DCM(2,1), et cetera
%
%
%   OUTPUTS
%       deltheta = profile of three orthogonal gyro outputs (body
%                  motion only) over time; for the i-th flight path
%                  segment, deltheta(i,1) = x-gyro output (nose positive),
%                  deltheta(i,2) = y-gyro output (right wing positive),
%                  deltheta(i,3) = z-gyro output (down positive)
%

%   REFERENCE:  Titterton, D. and J. Weston, STRAPDOWN
%               INERTIAL NAVIGATION TECHNOLOGY, Peter
%               Peregrinus Ltd. on behalf of the Institution
%               of Electrical Engineers, London, 1997, pp. 295-296.
%
%	M. & S. Braasch 01-98; modified 8-04 and 02-05
%	Copyright (c) 1998-2005 by GPSoft
%	All Rights Reserved.
%

if nargin<1,error('insufficient number of input arguments'),end

for k = 2:size(DCMnb_prof,1),
  dcmnb2=[DCMnb_prof(k,1:3); DCMnb_prof(k,4:6); DCMnb_prof(k,7:9)];
  dcmnb1=[DCMnb_prof(k-1,1:3); DCMnb_prof(k-1,4:6); DCMnb_prof(k-1,7:9)];
  A = dcmnb1 * dcmnb2';
  sigma_cross = logm(A);
  sigmaz(k-1,1) = real(sigma_cross(2,1));
  sigmay(k-1,1) = real(sigma_cross(1,3));
  sigmax(k-1,1) = real(sigma_cross(3,2));
end
deltheta = [sigmax sigmay sigmaz];

