function [Sat_XYZ,Pseudorange] = Seperate_Data(nSat,satsXYZ,prData,z)


%% Outputs
%%% Sat_XYZ -- Satetellites XYZ Coordinates in ECEF
%%% Pseudorange -- Pseudorange data 

%% Inputs
%%% nSat -- number of satellites in view
%%% satsXYZ --- XYZ of satellites in view in ECEF
%%% prData --- Pseudorange data 
%%% z --- loop index



%% Seperate Satellite Data 
Sat_XYZ(1:nSat(z),1:3,z) = satsXYZ(1:nSat(z),1:3,z);
%% Get PseudoRange Data 
Pseudorange(1:nSat(z),z) = prData(1:nSat(z),z);