function [Sat_XYZ_ENU,Unit_Vector_ENU] = Seperate_ENU_Data(nSat,satsXYZ,nomXYZ,Orgin,z)

%% Seperate Satellite Data 
    k=1;
    
    Sat_XYZ_ENU(1:nSat(z),1:3,z) = satsXYZ(1:nSat(z),1:3,z);
    
    for w = 1:nSat(z);
        k=(k-1)+1;
        Sat_XYZ_ENU(k,1:3,z) = xyz2enu(Sat_XYZ_ENU(k,1:3,z),Orgin);
    end

for j=1:nSat(z);
   
%     if j == 1; 
%         Unit_Vector_ENU(1,:,z)=(Sat_XYZ_ENU(1,:,z)-nomXYZ)/norm(Sat_XYZ_ENU(1,:,z)-nomXYZ);
%         EL(1,1:nSat(z)) = arcsin(dot(Unit_Vector_ENU(1,:,z)),(Sat_XYZ_ENU(1:nSat(z),3,z)));
%     else 
        y=j+1;
        Unit_Vector_ENU(y-1,:,z)=(Sat_XYZ_ENU(y-1,:,z)-nomXYZ)/norm(Sat_XYZ_ENU(y-1,:,z)-nomXYZ);
       % EL(y-1,1:nSat(z)) = arcsin(dot(Unit_Vector_ENU(y-1,1:3,z)),(Sat_XYZ_ENU((y-1),3,z)));
         AZ(y-1,1:nSat(z)) = 
    end 
end