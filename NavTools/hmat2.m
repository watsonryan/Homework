function h = hmat2(svmat,usrpos)
%HMAT2	Compute Direction Cosine Matrix
%
%	h = hmat2(svmat,usrpos)
%
%   INPUTS
%	svmat = matrix of satellite positions in
%	        user defined cartesian coordinates
%		   svmat(i,1:3) = x,y,z coordinates for satellite i
%	usrpos = estimated user position in user defined 
%	         cartesian coordinates
%
%   OUTPUTS
%	h = direction cosine matrix for GPS positioning
%
%   NOTE
%      This differs from HMAT.M in that this function will
%      output a matrix even for less than 4 satellites

%	Reference:
%                   Understanding GPS: Principles and Applications,
%	            Elliott D. Kaplan, Editor, Artech House Publishers,
%	            Boston, 1996.
%
%	Copyright (c) 1996-2005 by GPSoft
%
	[N,dum] = size(svmat);
	   tmppos = usrpos;
	   [m,n] = size(tmppos);
	   if m > n, tmppos = tmppos'; end,
	   h = ones(N,4);
	   for i = 1:N,
	       tmpvec = tmppos - svmat(i,:);
	       h(i,1:3) = tmpvec./norm(tmpvec);
	   end,

   	