%  sndemo30.m    Compute and plot autocorrelation of C/A PRN 10
%  
clear all
close all
%
[ca,g1,g2]=prncode(10);

numlag = 25;                      % To save time, we'll only
                                  % compute the first 25 lags
[lag,r] = gpscor(ca,ca,numlag);

plot(lag,r)
title('C/A PRN 10')
ylabel('autocorrelation value')
xlabel('lag in chips')
