function delta=tropgen(usrxyz,svxyz,humid)
%TROPGEN	Generate tropospheric delay in meters
%
%	delta = tropgen(usrxyz,svxyz,humid)
%
%   INPUTS
%	usrxyz(1:3) = true user position in ECEF cartesian coordinates.
%	svxyz(1:3) = position of satellite in ECEF cartesian coordinates.
%   humid = optional argument.  Humidity in percentage.
%
%   OUTPUTS
%	delta = tropospheric delay in meters corresponding to user and
%               satellite positions specified by usrxyz and svxyz

%	References: 
%                   Global Positioning System - Theory and Practice, 2nd ed.
%                   B. Hofmann-Wellenhof, H. Lichtenegger and J. Collins
%                   Springer-Verlag, Wien, New York, 1992.
%
%                   Understanding GPS: Principles and Applications,
%	            Elliott D. Kaplan, Editor, Artech House Publishers,
%	            Boston, 1996.
%
%                   Alfred Leick, GPS Satellite Surveying, 2nd ed.,
%	            Wiley-Interscience, John Wiley & Sons, 
%	            New York, 1995.
%
%	M. & S. Braasch 11-96
%	Copyright (c) 1996 by GPSoft
%	All Rights Reserved.
%

if nargin<2,error('insufficient number of input arguments'),end
if nargin<3,humid=50;end

% Atmospheric Parameters
p=1013; %std pressure in mbar
T=288.15; %temp in Kelvin
hp=0;
ht=1;
hh=0;
%
usrllh=xyz2llh(usrxyz);
hk=0.001*usrllh(3);
svenu=xyz2enu(svxyz,usrxyz);
beta=atan2(svenu(3),norm(svenu(1:2)));
%
ae=6378.137;
x=T-6.5*(hh-ht);
y=(7.5*(x-273.15))/(237.3+x-273.15);
eo=0.0611*humid*10^y;
tk=T+6.5*ht;
em=5.2459587;
es=eo*(tk/x)^(4*em);
te=T+6.5*(ht-hp);
psea=p*(tk/te)^em;
Nd0=77.624*psea/tk;
Nw0= (-12.92 + 371900/tk)*es/tk;

h(1)=1e-3*0.011385*psea/(Nd0*1e-6);
h(2)=1e-3*(0.011385/(Nw0*1e-6))*(1255/tk + 0.05)*es;
N1=Nd0*((h(1)-hk)/h(1))^4;
N2=Nw0*((h(2)-hk)/h(2))^4;

for i=1:2,
   a(i)=-(sin(beta)/(h(i)-hk));
   b(i)=-(1/(2*ae))*(((cos(beta))^2)/(h(i)-hk));
   alpha(1,i)=1;
   alpha(2,i)=4*a(i);
   alpha(3,i)=6*a(i)^2+4*b(i);
   alpha(4,i)=4*a(i)*(a(i)^2+3*b(i));
   alpha(5,i)=a(i)^4+12*a(i)^2*b(i)+6*b(i)^2;
   alpha(6,i)=4*a(i)*b(i)*(a(i)^2+3*b(i));
   alpha(7,i)=b(i)^2*(6*a(i)^2+4*b(i));
   alpha(8,i)=4*a(i)*b(i)^3;
   alpha(9,i)=b(i)^4;
   R(i)=sqrt( (ae+h(i))^2-(ae+hk)^2*(cos(beta))^2 )-((ae+hk)*sin(beta));
end,

sumr1=0;
sumr2=0;
for j=1:9,
   sum1=(alpha(j,1)/j)*R(1)^j;
   sumr1=sumr1+sum1;
   sum2=(alpha(j,2)/j)*R(2)^j;
   sumr2=sumr2+sum2;
end,
delta=1e-3*((N1*sumr1) + (N2*sumr2));

