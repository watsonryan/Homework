%% Ryan Watson
%% MAE 593
%% Homework #1
%% Part 1
%% Due 09/11/2014

clear all %% Clear data
clc   %% Clear workspace
load('dataSet1.mat') % Load Data
z=0; j=1; y=1; %% Set counters

Length = length(nSat);
c = 299792458; %% speed of light 

Sat_XYZ = zeros(max(nSat),3,Length); %% Allocate memory
Pseudorange = zeros(max(nSat),Length); %% Allocate memory
Computed_Pseudorange = zeros(max(nSat),1,Length); %% Allocate memory
Unit_Vector = zeros(max(nSat),3,Length); %% Allocate memory
XYZ_Estimate = zeros(3,Length); %% Allocate memory


for i = 1:Length-1;
    z=i+1; %% Set Counter
    
    %% Seperate Satellite Data 
    Sat_XYZ(1:nSat(z),1:3,z) = satsXYZ(1:nSat(z),1:3,z);
    %% Get PseudoRange Data 
    Pseudorange(1:nSat(z),z) = prData(1:nSat(z),z);


for j=1:nSat(z);
   
    if j == 1; 
        Computed_Pseudorange(1,1,z)=norm(Sat_XYZ(1,:,z)-nomXYZ)+clockBiasNom*c;
        Unit_Vector(1,:,z)=(Sat_XYZ(1,:,z)-nomXYZ)/norm(Sat_XYZ(1,:,z)-nomXYZ);
    else 
        y=j+1;
        Computed_Pseudorange(y-1,1,z)=norm(Sat_XYZ(y-1,:,z)-nomXYZ)+clockBiasNom*c;
        Unit_Vector(y-1,:,z)=(Sat_XYZ(y-1,:,z)-nomXYZ)/norm(Sat_XYZ(y-1,:,z)-nomXYZ);

    end 
end

G(1:nSat(z),:,z) = horzcat(Unit_Vector(1:nSat(z),:,z),ones(nSat(z),1));
deltaRho(1:nSat(z),:,z) = Computed_Pseudorange(1:nSat(z),:,z)-prData(1:nSat(z),z);
dX(:,:,z) = inv(G(:,:,z)'*G(:,:,z))*G(:,:,z)'*deltaRho(:,:,z);
XYZ_Estimate(:,z) = nomXYZ'+dX(1:3,:,z);
Error(:,z) = norm(XYZ_Estimate(:,z)-truthXYZ(:,z));
end

Initial_Error=norm(nomXYZ'-truthXYZ(:,1));
Final_Error= Error(Length);
Clock_Bias_Estimate = clockBiasNom+-1*dX(4)/c;
Initail_Clock_Error=norm(clockBiasNom-truthClockBias(1,1));
Final_Clock_Error=norm(Clock_Bias_Estimate-truthClockBias(1,Length));
solStr=sprintf('Initial 3D Error %.2f m, Final 3D Error %.2f m',Initial_Error,Final_Error);
disp(solStr)
% report clock estimation performance in meters
clockSolStr=sprintf('\nInitial Clk Error %.4f m, Final Clk Error %.4f m\n',...
    Initail_Clock_Error,Final_Clock_Error);
disp(clockSolStr)

