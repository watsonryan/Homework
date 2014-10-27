function svxyz = svpos(r,toe,Mo,OMGo,incl,t)
%SVPOS  Compute satellite position given Kepler parameters
%	corresponding to ideal circular orbits.  WGS-84
%
%	svxyz = svpos(r,toe,Mo,OMGo,incl,t)
%
%   INPUTS
%	r = orbit radius (semi-major axis) in meters
%	toe = reference time for Kepler parameters (time of ephemeris)
%	      in seconds
%	Mo = Mean anomaly at reference time in degrees
%	OMGo = Longitude of the ascending node at weekly epoch in degrees
%	incl = inclination angle of orbital plane in degrees
%	t = time at which to evaluate satellite position (in seconds)
%
%   OUTPUTS
%	svxyz(1,1) = ECEF x-coordinate of satellite in meters
%	svxyz(2,1) = ECEF y-coordinate of satellite in meters
%	svxyz(3,1) = ECEF z-coordinate of satellite in meters

%	Reference: Understanding GPS: Principles and Applications,
%	           Elliott D. Kaplan, Editor, Artech House Publishers,
%	           Boston, 1996.
%
%	M. & S. Braasch 10-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.

	mu = 3986005e8;
	OMGedot = 7.2921151467e-5;
	inclr = incl*pi/180;
	n = sqrt(mu/(r*r*r));
	tk = t - toe;
	Mk = Mo*pi/180 + n*tk;
	OMGk = OMGo*pi/180 - OMGedot*tk - OMGedot*toe;
	ip(1,1) = r*cos(Mk);
	ip(2,1) = r*sin(Mk);
	
	cosomg = cos(OMGk);
	sinomg = sin(OMGk);
	cosi = cos(inclr);
	ROT = [cosomg  -cosi*sinomg; ...
	       sinomg   cosi*cosomg; ...
	         0       sin(inclr)];
	svxyz = ROT*ip;
       	
