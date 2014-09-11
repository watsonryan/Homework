%% Ryan Watson
%% MAE 593
%% Homework #1
%% Part 2
%% Due 09/11/2014


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Inital Inputs  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all %% Clear data   
clc   %% Clear workspace
load('dataSet1.mat') % Load Data
z=1; j=1; y=1; %% Set counters
Orgin = nomXYZ; %% Orgin for converting to ENU
sigma_URE = sigma_URE^2*eye(4);

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
    [P,R,EL,W] = ELandAV(z,nSat,truthXYZ,Unit_Vector,W);
    [G,deltaRho,delta_x,delta_x_Weight,XYZ_Estimate,XYZ_Estimate_Weighted] = ECEF_Estimate(nSat,Unit_Vector, Computed_Pseudorange, prData,nomXYZ,z,W);
    H(:,:,z) = sigma_URE*(inv(G(:,:,z)'*G(:,:,z)));
    PDOP(:,z) = (sqrt(((H(1,1,z)^2)+(H(2,2,z)^2)+(H(3,3,z)^2))))/(sigma_URE(1,1));
    GDOP(:,z) = (sqrt(((H(1,1,z)^2)+(H(2,2,z)^2)+(H(3,3,z)^2) +(H(4,4,z)^2))))/((sigma_URE(1,1)));
    TDOP(:,z) = (sqrt((H(4,4,z)^2)))/(sigma_URE(1,1));  
    Error_ECEF(:,z) = norm(XYZ_Estimate(:,z)-truthXYZ(:,z));
    [Estimated_ENU(:,:,z),R_Estimate(:,:,z)] = xyz2enu(XYZ_Estimate(:,z)',Orgin);
    Error_ECEF_Weighted(:,z) = norm(XYZ_Estimate_Weighted(:,z)-truthXYZ(:,z));
    [True_ENU(:,:,z),R_True(:,:,z)] = xyz2enu(truthXYZ(:,z)',Orgin);
    Error_ENU(:,z) = norm(Estimated_ENU(:,:,z)-True_ENU(:,:,z));
    Clock_Bias_Estimate(:,z) = clockBiasNom-1*(delta_x(4)/Speed_of_Light);
    
    
    
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
    figure()
    plot(Error_ECEF-Error_ENU)
    title('Error_ECEF - Error_ENU')
    figure ()
    plot(PDOP,'x')
    title('PDOP')
    figure()
    plot(TDOP,'x')
    title('TDOP')
    figure()
    plot(GDOP,'x')
    title('GDOP')

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%  END PLOTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%