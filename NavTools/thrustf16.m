function thr = thrustf16(pow,alt,rmach)
%THRUSTF16    F-16 Engine thrust model
%
%  thr = thrustf16(pow,alt,rmach)
%
%  INPUTS
%     pow = power setting in percent
%     alt = altitude in feet
%     rmach = mach number
%
%  OUTPUT
%     thr = thrust
%
%  REFERENCE
%      The FORTRAN version of this program is given in:
%      AIRCRAFT CONTROL AND SIMULATION, 2nd edition, by B. Stevens
%      and F. Lewis, Wiley, Hoboken, NJ, 2003.

%	M. & S. Braasch 02-05
%	MATLAB version: Copyright (c) 2005 by GPSoft LLC
%	All Rights Reserved.
%
a = [1060    670   880   1140  1500 1860; ...
     635     425   690   1010  1330 1700; ...
     60       25   345    755  1130  1525; ...
     -1020  -710  -300    350   910  1360; ...
     -2700  -1900  -1300  -247  600  1100; ...
     -3600  -1400  -595  -342  -200 700];

% mil data
b = [12680  9150  6200  3950  2450  1400; ...
     12680  9150  6313  4040  2470  1400; ...
     12610  9312  6610  4290  2600   1560; ...
     12640  9839  7090  4660  2840  1660; ...
     12390  10176  7750  5320  3250  1930; ...
     11680  9848  8050  6100  3800  2310];

% max data
c = [20000  15000 10800 7000 4000 2500; ...
     21420 15700  11225  7323  4435  2600; ...
     22700  16860  12250  8154  5000  2835; ...
     24240  18910  13760  9285  5700  3215; ...
     26070  21075  15975  11115  6860  3950; ...
     28886  23319  18300  13484  8642  5057];

h = .0001*alt;
i = fix(h);
if i >= 5, i = 4; end
%%
if i <= 0, i = 0; end
%%
dh = h - i;
rm = 5*rmach;
m = fix(rm);
if m >= 5, m = 4; end
%%
if m <= 0, m = 0; end
%%
dm = rm - m;
cdh = 1 - dh;

i = i + 1;
m = m + 1;

s = b(m,i)*cdh + b(m,i+1)*dh;
t = b(m+1,i)*cdh + b(m+1,i+1)*dh;
tmil = s + (t-s)*dm;
if pow < 50,
    s = a(m,i)*cdh + a(m,i+1)*dh;
    t = a(m+1,i)*cdh + a(m+1,i+1)*dh;
    tidl = s + (t-s)*dm;
    thr = tidl + (tmil-tidl)*pow*0.02;
else
    s = c(m,i)*cdh + c(m,i+1)*dh;
    t = c(m+1,i)*cdh + c(m+1,i+1)*dh;
    tmax = s + (t-s)*dm;
    thr = tmil + (tmax-tmil)*(pow-50)*0.02;
end
    