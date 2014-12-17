%% Ryan Watson
%% MAE 593
%% Homework #1
%% Part 1
%% Constant Nominal Value
%% Due 09/11/2014

%clear all %% Clear data
clc   %% Clear workspace
load('dataSet3.mat') % Load Data
z=0; j=1; y=1; i=1;%% Set counters
Orgin = nomXYZ; %% Orgin for converting to ENU


T = 288.15;    % Kelvin
Tot_pressure = 1013;    %mbar
Par_pressure = 12.8;    %mbar



Length = length(nSat);
c = 299792458; %% speed of light 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  MEMORY ALLOCATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sat_XYZ = zeros(max(nSat),3,Length); 
Pseudorange = zeros(max(nSat),Length); 
Computed_Pseudorange = zeros(max(nSat),1,Length);   
Unit_Vector = zeros(max(nSat),3,Length); 
XYZ_Estimate = zeros(3,Length); 
nom_XYZ = zeros(3,Length); 
nom_XYZ(:,1) = nomXYZ';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1:Length;
    z=(i-1)+1; %% Set Counter
    
    if z ==1 
        
        %% Seperate Satellite Data 
        Sat_XYZ(1:nSat(1),:,1) = satsXYZ(1:nSat(1),:,1);
        %% Get PseudoRange Data 
        Pseudorange(1:nSat(1),1) = prDataP1(1:nSat(1),1);
    
    else
       
        %% Seperate Satellite Data 
        Sat_XYZ(1:nSat(z),:,z) = satsXYZ(1:nSat(z),:,z);
        %% Get PseudoRange Data 
        Pseudorange(1:nSat(z),z) = prDataP1(1:nSat(z),z);
    
    end
    
    %% Seperate Satellite Data 
    Sat_XYZ(1:nSat(z),:,z) = satsXYZ(1:nSat(z),:,z);
    %% Get PseudoRange Data 
    Pseudorange(1:nSat(z),z) = prDataP1(1:nSat(z),z);


for j=1:nSat(z);
   
    if j == 1; 
        Computed_Pseudorange(1,1,z)=norm(Sat_XYZ(1,:,z)-nomXYZ)+clockBiasNom*c;
        Unit_Vector(1,:,z)=(Sat_XYZ(1,:,z)-nomXYZ)/norm(Sat_XYZ(1,:,z)-nomXYZ);P(1,z) = sqrt(((truthXYZ(1,z))^2)+((truthXYZ(2,z))^2));
          
        
    else 
        y=j+1;
        Computed_Pseudorange(y-1,1,z)=norm(Sat_XYZ(y-1,:,z)-nomXYZ)+clockBiasNom*c;
        Unit_Vector(y-1,:,z)=(Sat_XYZ(y-1,:,z)-nomXYZ)/norm(Sat_XYZ(y-1,:,z)-nomXYZ);
         

    end 
end


G(1:nSat(z),:,z) = horzcat(Unit_Vector(1:nSat(z),:,z),ones(nSat(z),1));
deltaRho(1:nSat(z),:,z) = Computed_Pseudorange(1:nSat(z),:,z)-prDataP1(1:nSat(z),z);
dX(:,:,z) = inv(G(:,:,z)'*G(:,:,z))*G(:,:,z)'*deltaRho(:,:,z);
XYZ_Estimate(:,z) = nomXYZ'+dX(1:3,:,z);
nom_XYZ(:,z) = XYZ_Estimate(:,z);
Error_ECEF(:,z) = norm(XYZ_Estimate(:,z)-truthXYZ(:,z));
[Estimated_ENU(:,:,z),R_Estimate(:,:,z)] = xyz2enu(XYZ_Estimate(:,z)',Orgin);
[True_ENU(:,:,z),R_True(:,:,z)] = xyz2enu(truthXYZ(:,z)',Orgin);
Error_ENU(:,z) = norm(Estimated_ENU(:,:,z)-True_ENU(:,:,z));


llh(z,:) = xyz2llh(XYZ_Estimate(:,z)');

Td = .002277*(1+.0026*(cos(llh(z,1)))+(.0028*(llh(z,3))))*Tot_pressure;
Tw = .002277*((1255/T)+.05)*Par_pressure;
Tropo(:,z) = Td+Tw;

Tropo(:,z) = Td+Tw;
        
end

Initial_Error=norm(nomXYZ'-truthXYZ(:,1));
Final_Error= Error_ECEF(Length);
Clock_Bias_Estimate = clockBiasNom+-1*dX(4)/c;
Initail_Clock_Error=norm(clockBiasNom-truthClockBias(1,1));
Final_Clock_Error=norm(Clock_Bias_Estimate-truthClockBias(1,Length));
solStr=sprintf('Initial 3D Error %.2f m, Final 3D Error %.2f m',Initial_Error,Final_Error);
disp(solStr)
% report clock estimation performance in meters
clockSolStr=sprintf('\nInitial Clk Error %.4f m, Final Clk Error %.4f m\n',...
    Initail_Clock_Error,Final_Clock_Error);
disp(clockSolStr)
figure()
plot(Error_ECEF)
figure()
plot(Error_ENU)
figure()
plot(Error_ECEF-Error_ENU)

