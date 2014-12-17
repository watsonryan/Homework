%% Ryan Watson
%% MAE 593
%% Homework #1
%% Part 2
%% Due 09/11/2014


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Inital Inputs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all %% Clear data   
clc   %% Clear workspace
load('dataSet3.mat') % Load Data
z=1; j=1; y=1; %% Set counters
Orgin = nomXYZ; %% Orgin for converting to ENU
sigma_URE = sigma_URE^2*eye(4);


T = 288.15;    % Kelvin
Tot_pressure = 1013;    %mbar
Par_pressure = 12.8;    %mbar

prData = prDataP1;

Length = length(nSat);
Speed_of_Light = 299792458;  
[Sat_XYZ,Pseudorange,Computed_Pseudorange,Unit_Vector,XYZ_Estimate,PDOP,GDOP,TDOP,W,delta_x,delta_x_Weight] = Memory(nSat,Length);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  END INITIAL INPUTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:Length-1;
    
    z=(i-1)+1; 
    
    [Sat_XYZ,Pseudorange] = Seperate_Data(nSat,satsXYZ,prData,z);
    [Computed_Pseudorange,Unit_Vector] = Computed(Sat_XYZ,nomXYZ,clockBiasNom,nSat,Speed_of_Light,z);
    [P,R,EL,W,map] = ELandAV(z,nSat,truthXYZ,Unit_Vector,W);
    [G,deltaRho,delta_x,delta_x_Weight,XYZ_Estimate,XYZ_Estimate_Weighted] = ECEF_Estimate(nSat,Unit_Vector, Computed_Pseudorange, prData,nomXYZ,z,W);
    H(:,:,z) = sigma_URE*(inv(G(:,:,z)'*G(:,:,z))); 
    Error_ECEF(:,z) = norm(XYZ_Estimate(:,z)-truthXYZ(:,z));
    [Estimated_ENU(:,:,z),R_Estimate(:,:,z)] = xyz2enu(XYZ_Estimate(:,z)',Orgin);
    Error_ECEF_Weighted(:,z) = norm(XYZ_Estimate_Weighted(:,z)-truthXYZ(:,z));
    [True_ENU(:,:,z),R_True(:,:,z)] = xyz2enu(truthXYZ(:,z)',Orgin);
    Error_ENU(:,z) = norm(Estimated_ENU(:,:,z)-True_ENU(:,:,z));
    Clock_Bias_Estimate(:,z) = clockBiasNom-1*(delta_x(4)/Speed_of_Light);
    
    [Tropo] = tropo(z,Tot_pressure, Par_pressure,T, XYZ_Estimate);
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% END MAIN LOOP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%  PLOTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure()
    plot(Error_ECEF)
    title('ERROR_ECEF')
    figure()
    plot(Error_ENU)
    title('Error_ENU')
    figure()
    plot(Error_ECEF_Weighted)
    title('Weighted_ECEF')

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%  END PLOTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%