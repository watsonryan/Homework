%  sndemo29.m    Compute and plot autocorrelation of G1
%  
clear all
close all
%
[ca,g1,g2]=prncode(10);

numlag = 25;                      % To save time, we'll only
                                  % compute the first 25 lags
[lag,r] = gpscor(g1,g1,numlag);

plot(lag,r)
title('G1 code')
ylabel('autocorrelation value')
xlabel('lag in chips')
