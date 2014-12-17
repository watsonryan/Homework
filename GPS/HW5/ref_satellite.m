function [ref_prn,ref_sat_xyz] = ref_satellite(el_R,sd_prns,satsXYZ)


    ref_sat_candidate_vec=find(el_R==max(el_R));
    ref_sat_candidate_prn = sd_prns(ref_sat_candidate_vec,1);

    ref_sat_xyz(1,1:3,1) = [satsXYZ(1,ref_sat_candidate_prn,1)...
                           satsXYZ(1,ref_sat_candidate_prn,2)...
                           satsXYZ(1,ref_sat_candidate_prn,3)];
    ref_prn = ref_sat_candidate_prn;
    ref_prn_ = ref_sat_candidate_prn;
