function [clkcorr_m,grpdel] = svclkcorr(svid,towsec,Ek)
%SVCLKCORR	Compute satellite clock correction
%
%	[clkcorr_m,grpdel] = SVCLKCORR(svid,towsec,Ek)
%
%   INPUTS
%	svid = satellite identification number
%	towsec = GPS time of week in seconds
%   Ek = eccentric anomaly of the satellite in radians
%
%   OUTPUTS
%	clkcorr_m = satellite clock correction expressed in meters
%   grpdel = group delay of the satellite expressed in meters
%
%   NOTE: The group delay term, TGD, is NOT incorporated into
%         the clock correction term.  The single-frequency user
%         must implement this externally to this routine.  Dual
%         frequency users do not need to apply TGD.  See sections
%         20.3.3.3.3 through 20.3.3.3.3.3 of ICD-GPS-200 for 
%         more details

%	M. & S. Braasch     September 2002
%             Revised   February 2003
%	Copyright (c) 2002-2003 by GPSoft
%	All Rights Reserved.
%

global SV_ID_VEC TOC_YEAR TOC_MONTH TOC_DAY TOC_HOUR
global TOC_MINUTE TOC_SEC AF0 AF1 AF2
global TOE ECCEN SQRTSMA TGD

if nargin<2,error('insufficient number of input arguments'),end

toesv = TOE(svid);
if toesv == 604800,
    day_ambiguity = 0;
else
    day_ambiguity = floor(toesv/86400);
end

toc = TOC_HOUR(svid)*3600 + TOC_MINUTE(svid)*60 + TOC_SEC(svid) + day_ambiguity*86400;

tk = towsec - toc;
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

F = -4.442807633e-10;
deltatr = F*ECCEN(svid)*SQRTSMA(svid)*sin(Ek);

deltatsvL1 = AF0(svid) + AF1(svid)*tk + AF2(svid)*tk*tk + deltatr;
clkcorr_m = deltatsvL1*299792458;
grpdel = TGD(svid);