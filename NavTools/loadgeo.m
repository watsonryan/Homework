function loadgeo
%LOADGEO	GEOstationary satellite orbital parameters
%		Load the Kepler parameters for ideal circular orbits
%       from the matrix in geokep.mat and maintain as global variables
%
%	loadgeo
%
%   GLOBAL VARIABLES
%	SVIDV =	vector of satellite identification numbers
%	MV = 	vector of Mean anomalies for the satellites in svid
%		(at reference time) in degrees
%	RV =	vector orbit radii for the satellites in svid 
%		(semi-major axis) in meters
%	TOEV =	vector of reference times for the Kepler parameters 
%		for the satellites in svid (time of ephemeris) in seconds
%	OMGV = 	vector of longitudes of the ascending nodes for the 
%		satellites in svid (at weekly epoch) in degrees
%	INCLV = vector of inclination angles of orbital planes of
%		the satellites in svid (in degrees)

%	Copyright (c) 2003 by GPSoft
%
	global SVIDV MV OMGV RV INCLV TOEV

	load geokep

	SVIDV = [SVIDV; geokep(:,1)];
	MV = [MV; geokep(:,2)];
	OMGV = [OMGV; geokep(:,3)];
	RV = [RV; geokep(:,4)];
	INCLV = [INCLV; geokep(:,5)];
	TOEV = [TOEV; geokep(:,6)];  
     	
