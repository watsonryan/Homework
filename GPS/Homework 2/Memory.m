function [Sat_XYZ,Pseudorange,Computed_Pseudorange,Unit_Vector,XYZ_Estimate,PDOP,GDOP,TDOP,W,delta_x,delta_x_Weight] = Memory(nSat,Length)

%%% This function is used to allocate memory for certain variables

%% Inputs 
%%% nSat --- Number of satellites in view at epoch
%%% Length --- Number of epochs



Sat_XYZ = zeros(max(nSat),3,Length); %% Allocate memory
Pseudorange = zeros(max(nSat),Length); %% Allocate memory
Computed_Pseudorange = zeros(max(nSat),1,Length); %% Allocate memory
Unit_Vector = zeros(max(nSat),3,Length); %% Allocate memory
XYZ_Estimate = zeros(3,Length); %% Allocate memory
PDOP = zeros(1,Length); %% Allocate memory
GDOP = zeros(1,Length); %% Allocate memory
TDOP = zeros(1,Length); %% Allocate memory
W=zeros(max(nSat),max(nSat),Length); %% Allocate memory
delta_x=zeros(4,1,Length);
delta_x_Weight=zeros(4,1,Length);