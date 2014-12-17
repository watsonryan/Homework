function [G,deltaRho,delta_x,delta_x_Weight,XYZ_Estimate,XYZ_Estimate_Weighted] = ECEF_Estimate(nSat,Unit_Vector, Computed_Pseudorange, prData,nomXYZ,z,W)

%% Outputs
%%% G --- Geometry Matrix
%%% deltaRho -- 
%%% delta_x -- 
%%% XYZ_Estimate --- Estimated XYZ in ECEF


%% Inputs 
%%% nSat -- number of satellites in view
%%% Unit_Vector --- unit vector between user and satellite
%%% Computed_Pseudorange --- computed distance between user and satellite
%%% prData --- psuedorange data
%%% nomXYZ --- nominal XYZ in ECEF
%%% z --- loop index



G(1:nSat(z),:,z) = horzcat(Unit_Vector(1:nSat(z),:,z),ones(nSat(z),1));
deltaRho(1:nSat(z),:,z) = Computed_Pseudorange(1:nSat(z),:,z)-prData(1:nSat(z),z);
delta_x(:,:,z) = inv(G(:,:,z)'*W(1:nSat(z),1:nSat(z),z)*G(:,:,z))*G(:,:,z)'*W(1:nSat(z),1:nSat(z),z)*deltaRho(:,:,z);
delta_x_Weight(:,:,z) = inv(G(:,:,z)'*G(:,:,z))*G(:,:,z)'*deltaRho(:,:,z);
XYZ_Estimate(:,z) = nomXYZ'+delta_x(1:3,:,z);
XYZ_Estimate_Weighted(:,z) = nomXYZ'+delta_x_Weight(1:3,:,z);