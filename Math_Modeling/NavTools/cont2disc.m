function [Ad,Bd] = cont2disc(A,B,T)
%CONT2DISC		Discretization of state model. 
%       
%	[Ad,Bd] = CONT2DISC(A,B,T)
%
%   INPUTS
%       A = matrix relating state vector to its derivative
%       B = matrix relating input vector to the derivative
%           of the state vector
%       T = time-step in seconds
%
%   OUTPUTS
%       Ad, Bd = discrete-equivalents of A & B
%

%  REFERENCES
%       Kamen, E and B. Heck, FUNDAMENTALS OF SIGNALS AND SYSTEMS USING
%       MATLAB, Prentice-Hall, Upper Saddle River, NJ, 1997.
%
%	M. & S. Braasch    May 2004
%	Copyright (c) 2004 by GPSoft LLC
%	All Rights Reserved.
%

   if nargin<3,error('insufficient number of input arguments'),end
   
Ad = expm(A*T);

Bd = zeros(size(B));
for k = 0:10,
    Bd = Bd + (1/factorial(k+1))*(A^k)*(T^(k+1))*B;
end
