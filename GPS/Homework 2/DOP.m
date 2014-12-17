function [H,PDOP,GDOP,TDOP] = DOP(sigma_URE,G,z)

%% This function is used to calculate:
%% PDOP -- Position Dilution of Precision
%% GDOP -- Geometric Dilution of Precision
%% TDOP -- Time Dilution of Precision
%% H -- Observaion Matrix 

%% Inputs 
%%% Sigma_URE -- Variance
%%% G -- Geometry Matrix
%%% z -- loop index


H(:,:,z) = sigma_URE*(inv(G(:,:,z)'*G(:,:,z)));

PDOP(:,z) = sqrt(((H(1,1,z)^2)+(H(2,2,z)^2)+(H(3,3,z)^2)));
GDOP(:,z) = sqrt(((H(1,1,z)^2)+(H(2,2,z)^2)+(H(3,3,z)^2) +(H(4,4,z)^2)));
TDOP(:,z) = sqrt((H(4,4,z)^2));

