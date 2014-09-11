function [Computed_Pseudorange,Unit_Vector] = Computed(Sat_XYZ,nomXYZ,clockBiasNom,nSat,c,z);

%% Outputs
%%% Computed_Pseudorange -- Computed Pseudorange based upon raw data
%%% Unit_Vector -- Unit Vector between user and satellite


%% Inputs
%%% Sat_XYZ -- Satellite cooridnated in ECEF
%%% nomXYZ -- nominal value for linerization of model
%%% nSat -- number of satellites in view
%%% clockBiasNom -- nominal clock bias
%%% c -- speed of light (m/s)
%%% z -- loop index


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