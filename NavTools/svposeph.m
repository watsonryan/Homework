function [svxyz,Ek] = svposeph(svid,t)
%SVPOSEPH  Compute satellite position given the
%  broadcast ephemeris parameters.  WGS-84
%
%	[svxyz,Ek] = svposeph(svid,t)
%
%   INPUTS
%   svid = identification number of satellite
%	t = time at which to evaluate satellite position (in seconds)
%
%   OUTPUTS
%	svxyz(1,1) = ECEF x-coordinate of satellite in meters
%	svxyz(2,1) = ECEF y-coordinate of satellite in meters
%	svxyz(3,1) = ECEF z-coordinate of satellite in meters
%   Ek = Eccentric anomaly
%
%   NOTE:  SVPOSEPH assumes the user has used the function LOADRINEXN to
%   load the ephemeris parameters into global memory

%	Reference: ICD-GPS-200 ...
%
%	M. & S. Braasch  -  September 2002
%	Copyright (c) 2002 by GPSoft LLC
%	All Rights Reserved.

global IODE CRS DELTAN MZERO CUC ECCEN CUS SQRTSMA
global TOE CIC OMEGAZERO CIS IZERO CRC ARGPERI OMEGADOT
global IDOT TOE_WN

	mu = 3986005e8;
	OMGedot = 7.2921151467e-5;
    A = ( SQRTSMA(svid) )^2;
	n_o = sqrt(mu/(A*A*A));
    n = n_o + DELTAN(svid);
    tk = t - TOE(svid);
      iter=0;
      while tk > 302400,
         tk = tk - 604800;
         iter=iter+1;
         if iter > 3, error('Input time should be time of week in seconds'), end
      end
      iter=0;
      while tk < -302400,
         tk = tk + 604800;
         iter=iter+1;
         if iter > 3, error('Input time should be time of week in seconds'), end
      end
   Mk = MZERO(svid) + n*tk;
   
   Ek = Mk; sep = 1; oldEk = Ek;
   iter = 0;
   while sep > 1e-13,
      Ek = Mk + ECCEN(svid)*sin(Ek);
      sep = abs(Ek - oldEk);
      oldEk = Ek;
      iter = iter + 1;
      if iter > 10, break, end
   end
   
   sin_nu_k = ( (sqrt(1-ECCEN(svid)*ECCEN(svid)))*sin(Ek) )/( 1 - ECCEN(svid)*cos(Ek) );
   cos_nu_k = ( cos(Ek) - ECCEN(svid) )/( 1 - ECCEN(svid)*cos(Ek) );
   nu_k = atan2(sin_nu_k,cos_nu_k);
   
   PHIk = nu_k + ARGPERI(svid);
   
   c2 = cos(2*PHIk);
   s2 = sin(2*PHIk);
   delta_uk = CUS(svid)*s2 + CUC(svid)*c2;
   delta_rk = CRS(svid)*s2 + CRC(svid)*c2;
   delta_ik = CIS(svid)*s2 + CIC(svid)*c2;
   
   uk = PHIk + delta_uk;
   rk = A*(1-ECCEN(svid)*cos(Ek)) + delta_rk;
   ik = IZERO(svid) + delta_ik + IDOT(svid)*tk;
   
   ip(1,1) = rk*cos(uk);    % satellite position In orbital Plane
   ip(2,1) = rk*sin(uk);   
   
	OMGk = OMEGAZERO(svid) + (OMEGADOT(svid)-OMGedot)*tk - OMGedot*TOE(svid);
	
	cosomg = cos(OMGk);
	sinomg = sin(OMGk);
	cosi = cos(ik);
	ROT = [cosomg  -cosi*sinomg; ...
	       sinomg   cosi*cosomg; ...
	         0       sin(ik)];
	svxyz = ROT*ip;
       	
