% Andrew Saiko - MAE 593 - GPS - Homework #1
format long g
clear
clc
load('F:\MAE 593\HW1\dataSet1.mat')
% This loads the data provided for the homework
i=1;
j=1;
N=2400;
% Sets a counter for each iteration

%% Truth Inputs for each iteration, i
ClockBias_Truth = truthClockBias(i);
XYZ_Truth = truthXYZ(:,i);
%% Given Inputs
nomXYZ;
clockBiasNom;
c = 299792458; % m/s 
%% Organizing satellite Data
satellites_XYZ = satsXYZ(1:nSat(i),1:3,i);
% This grabs the XYZ information for each satellite for the number of 
% satellites defined
pseudoR_Data = prData(1:nSat(i),i);
% This grabs the pseudorange data for each satellite for the # of
% satellites defined
prComputed=zeros(nSat(i),1,N);
uNom2Sat=zeros(nSat(i),3,N);
% pre allocate sizes
%% Calculated the pseudorange and the unit vector to satellite location
for j=1:nSat(i);
    % calculate computed or expected psuedorange
    prComputed(j)=norm(satellites_XYZ(j,:)-nomXYZ)+clockBiasNom*c;
    % calculate the unit vectors from nominal position
    % to satellite location 
    % which make up the partials 
    % of the observation model
    uNom2Sat(j,:,i)=(satellites_XYZ(j,:)-nomXYZ)/norm(satellites_XYZ(j,:)-nomXYZ);
end
%%
% add the partial of the clock bias, 
% observation model
G = horzcat(uNom2Sat(:,:,i),ones(nSat(i),1));
% form the innovation vector
deltaRho = prComputed(:,:,i)-prData(1:nSat(i),i);
% find the delta solution
dX = inv(G'*G)*G'*deltaRho;
XYZ_Estimate = nomXYZ'+dX(1:3)
% the clock bias in has a negative sign that must be handled
clockBiasEst = clockBiasNom+-1*dX(4)
%% Determine Performance
initErr=norm(nomXYZ'-XYZ_Truth);
finalErr=norm(XYZ_Estimate-XYZ_Truth);

I_Clock_Error=norm(clockBiasNom-ClockBias_Truth);
F_Clock_Error=norm(clockBiasEst-ClockBias_Truth);
solStr=sprintf('Initial 3D Error %.2f m, Final 3D Error %.2f m',initErr,finalErr);
disp(solStr)
% report clock estimation performance in meters
clockSolStr=sprintf('\nInitial Clk Error %.4f m, Final Clk Error %.4f m\n',...
    I_Clock_Error,F_Clock_Error);
disp(clockSolStr)

