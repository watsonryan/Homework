function svxyz = svposalm(r,toe,Mo,OMGo,incl,t,eccen,argperi,omgdot,weeknum,currentweek)
%SVPOSALM  Compute satellite position given the
%  broadcast almanac parameters.  WGS-84
%
%	svxyz = svposalm(r,toe,Mo,OMGo,incl,t,eccen,argperi,omgdot,weeknum,currentweek)
%
%   INPUTS
%	r = orbit radius (semi-major axis) in meters
%	toe = reference time for Kepler parameters (time of ephemeris)
%	      in seconds
%	Mo = Mean anomaly at reference time in degrees
%	OMGo = Longitude of the ascending node at weekly epoch in degrees
%	incl = inclination angle of orbital plane in degrees
%	t = GPS time of week at which to evaluate satellite position (in seconds)
%  eccen = eccentricity of the orbit
%  argperi = argument of perigee (radians)
%  omgdot = rate of right ascension (radians/sec)
%
%   OUTPUTS
%	svxyz(1,1) = ECEF x-coordinate of satellite in meters
%	svxyz(2,1) = ECEF y-coordinate of satellite in meters
%	svxyz(3,1) = ECEF z-coordinate of satellite in meters
%
%   Note: If the almanac is given in the same week as the time
%         at which the satellite position is to be computed ('t'),
%         the almanac time of applicability must be less than
%         or equal to 't'

%	Reference: Understanding GPS: Principles and Applications,
%	           Elliott D. Kaplan, Editor, Artech House Publishers,
%	           Boston, 1996.
%
%	M. & S. Braasch 10-99; revised 10-02 and 08-03
%	Copyright (c) 1999-2003 by GPSoft LLC
%	All Rights Reserved.

if nargin < 10, weeknum=-999; end
if nargin < 11, currentweek=weeknum; end

	mu = 3986005e8;
	OMGedot = 7.2921151467e-5;
	inclr = incl*pi/180;
	n = sqrt(mu/(r*r*r));
    
    if currentweek == weeknum,
        tk = t - toe;
    else
        tk = t + (604800-toe) + 604800*(currentweek-weeknum-1);
    end
    
%     iter=0;
%     while tk > 302400,
%     while tk > 604800,
%         tk = tk - 604800;
%         iter=iter+1;
%         if iter > 3, error('Input time should be time of week in seconds'), end
%     end
%     iter=0;
%     while tk < -302400,
%     while tk < 0,
%         tk = tk + 604800;
%         iter=iter+1;
%         if iter > 3, error('Input time should be time of week in seconds'), end
%     end

     Mk = Mo*pi/180 + n*tk;
   
   Ek = Mk; sep = 1; oldEk = Ek;
   iter = 0;
   while sep > 1e-13,
      Ek = Mk + eccen*sin(Ek);
      sep = abs(Ek - oldEk);
      oldEk = Ek;
      iter = iter + 1;
      if iter > 10, break, end
   end
   
   sin_nu_k = ( (sqrt(1-eccen*eccen))*sin(Ek) )/( 1 - eccen*cos(Ek) );
   cos_nu_k = ( cos(Ek) - eccen )/( 1 - eccen*cos(Ek) );
   nu_k = atan2(sin_nu_k,cos_nu_k);
   
   PHIk = nu_k + argperi;
   rc = r*(1-eccen*cos(Ek));
	ip(1,1) = rc*cos(PHIk);
	ip(2,1) = rc*sin(PHIk);   
   
	OMGk = OMGo*pi/180 + (omgdot-OMGedot)*tk - OMGedot*toe;
	
	cosomg = cos(OMGk);
	sinomg = sin(OMGk);
	cosi = cos(inclr);
	ROT = [cosomg  -cosi*sinomg; ...
	       sinomg   cosi*cosomg; ...
	         0       sin(inclr)];
	svxyz = ROT*ip;
       	
