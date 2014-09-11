function [P,R,EL,W] = ELandAV(z,nSat,truthXYZ,Unit_Vector,W)
        
        P(1,z) = sqrt(((truthXYZ(1,z))^2)+((truthXYZ(2,z))^2));
        R(1,z) = sqrt(((truthXYZ(1,z))^2)+((truthXYZ(2,z))^2)+((truthXYZ(3,z))^2));
    
        East(z,:) = [(-truthXYZ(2,z)^2)/P(1,z) (truthXYZ(1,z)^2)/P(1,z) 0 ];
        North(z,:) = [(((-truthXYZ(1,z))*(-truthXYZ(3,z)))/(P(1,z)*R(1,z)))  (((-truthXYZ(2,z))*(-truthXYZ(3,z)))/(P(1,z)*R(1,z)))  (P(1,z)/R(1,z))];
        Up(z,:) = [((truthXYZ(1,z))/(R(1,z)))  (truthXYZ(2,z)/(R(1,z)))  (truthXYZ(3,z)/R(1,z))];

        Up_Ones = ones(nSat(z),3);
        Up_z(:,:,z)=[Up(z,1)*Up_Ones(:,1)  Up(z,2)*Up_Ones(:,2)  Up(z,3)*Up_Ones(:,3)];
        
        k=1;
        for k=1:nSat(z)
            
            EL(:,:,z) = dot(Unit_Vector(:,:,z),Up_z(:,:,z),2);
            W(k,k,z)=1/(diag(EL(k,:,z)));
            
            k=k+1;
        end
        
end
   


    
    

