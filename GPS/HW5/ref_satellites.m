function [ref_sat_prn] = ref_satellites(prn_prev,prn_current,...
                                        ref_prn_,ref_prn,satsXYZ,i,...
                                         xyz_wvmo);


if ~isempty(find(prn_current == ref_prn_, 1))   %% is old satellite still in view  
        
        vec_position = find(prn_current == ref_prn_); 

        xyz = xyz2enu([satsXYZ(i,ref_prn_,1),...
                        satsXYZ(i,ref_prn_,2),...
                        satsXYZ(i,ref_prn_,3)],...
                        xyz_wvmo);

        New_el=atan2(xyz(3),sqrt(xyz(1)^2+xyz(2)^2))*(180*pi);
    
        if New_el > 12.5     %% If in view and has elevation angle above 20 deg

            ref_sat_prn = ref_prn_;

        else                %% If in view but has elevation angle below 20 deg
            
            ref_sat_prn = ref_prn;
            
        end
    
else     %% If sat not still in view 
        
     ref_sat_prn = ref_prn;
        
end  
    
    
    clear temp