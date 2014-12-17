function tropodel = tropocorr(svxyz,usrxyz)
%TROPOCORR	Compute tropospheric correction for
%           the dry component using the Hopfield model
%
%	tropodel = tropocorr(svxyz,usrxyz)
%
%   INPUTS
%	svxyz = satellite position expressed in ECEF cartesian coordinates
%	usrxyz = user position expressed in ECEF cartesian coordinates
%
%   OUTPUTS
%	tropodel = estimate of tropospheric error (meters)



svenu = xyz2enu(svxyz,usrxyz);
el = atan2(svenu(3),norm(svenu(1:2)));

usrllh = xyz2llh(usrxyz);
h = usrllh(3)/1000;

%tropodel = 2.47/(sin(el)+0.0121);
tropodel = 2.4224*exp(-0.13345*h)/(sin(el)+0.026);