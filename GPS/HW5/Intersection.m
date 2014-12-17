function   [prn,nSD] = Intersection(R1_wv,R2_wv,R1_pa, R2_pa,...
                                         P1_wv,P2_wv,P1_pa,P2_pa,i)

    %Range Data 
    wv_1 = isnan(R1_wv); wv_1 = ~wv_1; SatID_R1_wv = find(wv_1 == 1);
    wv_2 = isnan(R2_wv); wv_2 = ~wv_2; SatID_R2_wv = find(wv_2 == 1);  
    pa_1 = isnan(R1_pa); pa_1 = ~pa_1; SatID_R1_pa = find(pa_1 == 1); 
    pa_2 = isnan(R2_pa); pa_2 = ~pa_2; SatID_R2_pa = find(pa_2 == 1); 
    
    %Phase Data
    wv_1 = isnan(P1_wv); wv_1 = ~wv_1; SatID_P1_wv = find(wv_1 == 1);  
    wv_2 = isnan(P2_wv); wv_2 = ~wv_2; SatID_P2_wv = find(wv_2 == 1);  
    pa_1 = isnan(P1_pa); pa_1 = ~pa_1; SatID_P1_pa = find(pa_1 == 1);  
    pa_2 = isnan(P2_pa); pa_2 = ~pa_2; SatID_P2_pa = find(pa_2 == 1); 
    
    id1=intersect(SatID_R1_wv,SatID_R1_pa); id2 = intersect(SatID_R2_pa,SatID_R2_wv); 
    id3 = intersect(id1,id2);
    id4=intersect(SatID_P1_wv,SatID_P1_pa); id5 = intersect(SatID_P2_pa,SatID_P2_wv); 
    id6 = intersect(id4,id5);
    
    prn = intersect(id3,id6);

    nSD=length(prn);