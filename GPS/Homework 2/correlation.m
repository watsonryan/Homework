function [C,Max] = correlation()


clear all
clc


load('mysteryCA(1).mat');
C = zeros(37,(2*1014)-1);
PRN=zeros(37,1023);
PRN=cacode(1:37);
Max=zeros(37,1);


y=0;
k=1;

for i=1:37
    
    y=y+1;
    C(y,:)= xcorr(mysteryCA(1:1014),PRN(y,10:1023),'coeff');
     
end
end