function [ca,g1,g2]=prncode(prnnum)
%PRNCODE		Generate C/A-code sequence
%
%	[CA,G1,G2] = PRNCODE(PRNNUM)
%
% INPUTS
%     prnnum = GPS PRN signal number (must be an integer in the range 1:37)
%
% OUTPUTS
%   ca = C/A-code sequence realized with +1's and -1's (digital 0's and 1's).
%        CA is the product of G1 and G2 (after G2 has been shifted)
%     g1, g2 = maximal length sequences with g2 shifted according to the
%              specification in GPS-ICD-200.

%	References: 
%                   Global Positioning System Standard Positioning Service
%                   Signal Specification; 2nd edition, June 5, 1995.
%
%                   ICD-GPS-200, NAVSTAR GPS Space Segment/Navigation User
%                   Interfaces (Public Release Version), ARINC Research
%                   Corporation, 11770 Warner Ave., Suite 210, Foutain
%                   Valley, CA 92708, July 3, 1991.
%
%                   Understanding GPS: Principles and Applications,
%	            Elliott D. Kaplan, Editor, Artech House Publishers,
%	            Boston, 1996.
%

%	M. & S. Braasch 12-96;  Revised 10-99;  Revised 07-02
%	Copyright (c) 1996 - 2002 by GPSoft LLC
%	All Rights Reserved.
%

if nargin<1,error('insufficient number of input arguments'),end
if ((prnnum<0)|(prnnum>37)),error('invalid PRN number.  Must be in the range 1:37'),end
if round(prnnum)~=prnnum,error('invalide PRN number.  Must be an integer'),end

codedelay = [5 6 7 8 17 18 139 140 141 251 252 254 255 256 257 258 469 470 471 472 473 474 ...
      509 512 513 514 515 516 859 860 861 862 863 950 947 948 950];

R1 = -1*ones(1,10);
R2 = R1;

for k=1:1023,
   G1(k)=R1(10);
   temp1 =R1(3)*R1(10);
   R1(2:10)=R1(1:9);
   R1(1)=temp1;
   
   G2(k)=R2(10);
   temp2 =R2(2)*R2(3)*R2(6)*R2(8)*R2(9)*R2(10);
   R2(2:10)=R2(1:9);
   R2(1)=temp2;   
end,

shift=codedelay(prnnum);
shiftG2(shift+1:1023)=G2(1:1023-shift);
shiftG2(1:shift)=G2(1023-shift+1:1023);

ca=G1.*shiftG2;
g1 = G1;
g2 = shiftG2;
