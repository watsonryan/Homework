%  sndemo28.m    Compute and plot part of C/A PRN 10
%  
clear all
close all
%
[ca,g1,g2]=prncode(10);

subplot(311)
plot(g1(1:20))
ylabel('G1')
axis([0 20 -2 2])

subplot(312)
plot(g2(1:20))
ylabel('G2')
axis([0 20 -2 2])

subplot(313)
plot(ca(1:20))
ylabel('C/A')
axis([0 20 -2 2])
