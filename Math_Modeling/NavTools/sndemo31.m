%  sndemo31.m    Compute and plot cross-correlation 
%                between C/A PRN 10 and C/A PRN 20
%  
clear all
close all
%
ca10=prncode(10);

ca20 = prncode(20);

numlag = 25;                      % To save time, we'll only
                                  % compute the first 25 lags
[lag,r] = gpscor(ca10,ca20,numlag);

plot(lag,r)
title('C/A PRN 10 and PRN 20')
ylabel('cross-correlation value')
xlabel('lag in chips')
axis([0 25 -200 1200])