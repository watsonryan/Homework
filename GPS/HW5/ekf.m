clear all; close all; clc; 

load('realDataMat')

epochs=288;
R1_wv = wvmo.C1;
R2_wv = wvmo.P2;
R1_pa = pagw.C1;
R2_pa = pagw.P2;
lambdaL1=299792458/1575.42e6;
lambdaL2=299792458/1227.6e6;
lambda1=299792458/1575.42e6;
lambda2=299792458/1227.6e6;
P1_wv = wvmo.PHASE1*299792458/1575.42e6;
P2_wv = wvmo.PHASE2*299792458/1227.6e6;
P1_pa = pagw.PHASE1*299792458/1575.42e6;
P2_pa = pagw.PHASE2*299792458/1227.6e6;
Sats_Used = zeros(12,288);
Used_Sats_XYZ = zeros(22,3, 288);


xyz_wvmo=llh2xyz([39.64589002*pi/180,-79.96992114*pi/180,326.1251]);
xyz_pagm=llh2xyz([39.89893047*pi/180,-80.16011165*pi/180,262.7526]);

xyz_baseline=xyz_pagm;
enu_baseline = xyz2enu(xyz_baseline,xyz_wvmo)';


idSVdd=zeros(12,288);

for i = 1:288%length(satTOWSEC)
    %Range Data 
    wv_1(:,i) = isnan(R1_wv(:,i)); wv_1(:,i) = ~wv_1(:,i); SatID_R1_wv = find(wv_1(:,i) == 1); numSats_R1_wv(i) = length(SatID_R1_wv);  
    wv_2(:,i) = isnan(R2_wv(:,i)); wv_2(:,i) = ~wv_2(:,i); SatID_R2_wv = find(wv_2(:,i) == 1); numSats_R2_wv(i) = length(SatID_R2_wv); 
    pa_1(:,i) = isnan(R1_pa(:,i)); pa_1(:,i) = ~pa_1(:,i); SatID_R1_pa = find(pa_1(:,i) == 1); numSats_R1_pa(i) = length(SatID_R1_pa); 
    pa_2(:,i) = isnan(R2_pa(:,i)); pa_2(:,i) = ~pa_2(:,i); SatID_R2_pa = find(pa_2(:,i) == 1); numSats_R2_pa(i) = length(SatID_R2_pa);
    
    %Phase Data
    wv_1(:,i) = isnan(P1_wv(:,i)); wv_1(:,i) = ~wv_1(:,i); SatID_P1_wv = find(wv_1(:,i) == 1); numSats_P1_wv(i) = length(SatID_P1_wv); 
    wv_2(:,i) = isnan(P2_wv(:,i)); wv_2(:,i) = ~wv_2(:,i); SatID_P2_wv = find(wv_2(:,i) == 1); numSats_P2_wv(i) = length(SatID_P2_wv); 
    pa_1(:,i) = isnan(P1_pa(:,i)); pa_1(:,i) = ~pa_1(:,i); SatID_P1_pa = find(pa_1(:,i) == 1); numSats_P1_pa(i) = length(SatID_P1_pa); 
    pa_2(:,i) = isnan(P2_pa(:,i)); pa_2(:,i) = ~pa_2(:,i); SatID_P2_pa = find(pa_2(:,i) == 1); numSats_P2_pa(i) = length(SatID_P2_pa);
    
    id1=intersect(SatID_R1_wv,SatID_R1_pa); id2 = intersect(SatID_R2_pa,SatID_R2_wv); id3 = intersect(id1,id2);
    id4=intersect(SatID_P1_wv,SatID_P1_pa); id5 = intersect(SatID_P2_pa,SatID_P2_wv); id6 = intersect(id4,id5);
    
    id = intersect(id3,id6);
    
    Number_of_Sats_Used(i) = length(id);
    Sats_Used(1:length(id),i) = id;

    nSD=length(id);
    
    for k=1:nSD
        
        ID_SD=id(k);
        
        sd_rangeL1(k,i) = -R1_wv(ID_SD,i)+R1_pa(ID_SD,i);
        sd_phaseL1(k,i) = -P1_wv(ID_SD,i)+P1_pa(ID_SD,i);
        sd_rangeL2(k,i) = -R2_wv(ID_SD,i)+R2_pa(ID_SD,i);
        sd_phaseL2(k,i) = -P2_wv(ID_SD,i)+P2_pa(ID_SD,i);
        tempENU=xyz2enu([satsXYZ(i,ID_SD,1),satsXYZ(i,ID_SD,2),satsXYZ(i,ID_SD,3)], xyz_wvmo);
        el_R(k)=atan2(tempENU(3),sqrt(tempENU(1)^2+tempENU(2)^2));%*(180/pi);
        
    end
   
    Angle = max(el_R)*(180/pi);
    ref_sat=find(el_R==max(el_R));
    
    if i == 1 
        
        id_ref=Sats_Used(idx(1),1);
        RefSat_XYZ(1,1:3,1) = [satsXYZ(1,id_ref,1) satsXYZ(1,id_ref,2) satsXYZ(1,id_ref,3)];
        id_ref_ = id_ref;
    end
    
    if ~isempty(find(Sats_Used(:,i) == id_ref_))   %% Still available 
        
        xxxxx = find(Sats_Used(:,i) == id_ref_); 
        newTemp=xyz2enu([satsXYZ(i,Sats_Used(xxxxx,i),1),satsXYZ(i,Sats_Used(xxxxx,i),2),satsXYZ(i,Sats_Used(xxxxx,i),3)], xyz_wvmo);
        New_el=atan2(newTemp(3),sqrt(newTemp(1)^2+newTemp(2)^2))*(180*pi);
    
    if i > 1 && New_el > 20
        ref_sat = find(Sats_Used(:,i) == idr(i-1)); 
    else 
        ref_sat = ref_sat;
    end
    
    else 
     ref_sat = ref_sat;
    end  
    
    temp=[Sats_Used(1:(ref_sat-1),i);...
        Sats_Used((ref_sat+1):Number_of_Sats_Used(i),i)];
    idSVdd(1:length(temp),i)=temp;
    
    clear temp
    
    z=0;
    for j=1:nSD
        ID_DD=id(j);
        if(ID_DD~=id_ref)
            dd_range1(j-z,i)=sd_rangeL1(j,i)-sd_rangeL1(idx(i));
            dd_range2(j-z,i)=sd_rangeL2(j,i)-sd_rangeL2(idx(i));
            dd_phase1(j-z,i)=sd_phaseL1(j,i)-sd_phaseL1(idx(i));
            dd_phase2(j-z,i)=sd_phaseL2(j,i)-sd_phaseL2(idx(i));
            
            Used_Sats_XYZ(j-z,1:3,i) = [satsXYZ(i,ID_DD,1) satsXYZ(i,ID_DD,2) satsXYZ(i,ID_DD,3)];
        else
            z=1;
        end
        
    end
    
    NumDD(i) = length(id)-1;
    Used_Sats_XYZ(NumDD+1:NumDD+NumDD,1:3,i)=Used_Sats_XYZ(1:NumDD,1:3,i); 
    el_R_ = el_R;
    el_R =[];
    ref_sat_=ref_sat;
    id_ref_ = id_ref;
    idx(i) = ref_sat;
    idr(i) = id_ref;
    
clear id
end

%% good
tru_usrxyzO=xyz_wvmo;
reference_satellite=idx(1);
nParm=3+2*NumDD(1);
nNonBiasParm=3;

%% good
x=[xyz_pagm-xyz_wvmo, dd_range1(1:NumDD(1),1)'/lambda1,dd_range2(1:NumDD(1),1)'/lambda2]';
%x=[0,0,0, dd_range1(1:NumDD(11),11)',dd_range2(1:NumDD(11),11)']';
x_fixed=zeros(nParm,epochs);

notFixed=0;
Sigma=3e3;

Px=diag([1,1,1,ones(1,NumDD(1)*2)*1000^2]);
%Px=diag([1,1,1,ones(1,NumDD(11)*2)*15^2]);

rval=0.05;
qval=0.1;
R=(rval^2)*(ones(NumDD(1)*2,NumDD(1)*2)+eye(NumDD(1)*2));
%R=(.02^2)*(ones(NumDD(11)*2,NumDD(11)*2)+eye(NumDD(11)*2));

%% 
x_save=zeros(200,epochs);
x_=x;
Px_=Px;
nParm_=nParm;
for i=2:epochs
    
    nDD = length(nonzeros(dd_range1(:,i)));
    reference_satellite=idx(i);
    nObs=nDD*2;
    nParm=3+2*nDD;
    nNonBiasParm=3;
    
    if i >2 && reference_satellite_ ~= reference_satellite;
        Amb_L1 = nonzeros(dd_range1(:,i))./lambdaL1;
        Amb_L2 = nonzeros(dd_range2(:,i))./lambdaL2;
        x = [x_(1:3); Amb_L1; Amb_L2];
        
        Px = Sigma*eye(2*nDD+3);
        Px(1:3,1:3) = Px_(1:3,1:3);
        
    else
        [x,Px]=addRemoveBiasParm(x_(1:nParm_),Px_,nNonBiasParm,nonzeros(idSVdd(:,i-1)),nonzeros(idSVdd(:,i)), Sigma^2);
    end
    
    
    R=(rval^2)*(ones(nDD*2,nDD*2)+eye(nDD*2));
    Q=zeros(nParm);
    for j=1:3
        Q(j,j)=qval;
    end
    
    Rsat=[satsXYZ(i,Sats_Used(idx(i),i),1), satsXYZ(i,Sats_Used(idx(i),i),2), satsXYZ(i,Sats_Used(idx(i),i),3)];
    vec2RSat=((Rsat-tru_usrxyzO)/norm(Rsat-tru_usrxyzO));
    uSats=[];
    for j=1:nDD
        xyzSat=[satsXYZ(i,idSVdd(j,i),1),satsXYZ(i,idSVdd(j,i),2),satsXYZ(i,idSVdd(j,i),3)];
        vecXYZSat=((xyzSat-tru_usrxyzO)/norm(xyzSat-tru_usrxyzO));
        uSats(j,1:3)=(vecXYZSat-vec2RSat);
    end
    uSats=[uSats;uSats]; 

    H=zeros(2*nDD,3+2*nDD);
    H(:,1:3)=uSats(:,1:3);
    for j=1:2*nDD
        if(j<=nDD)
           H(j,j+3)=lambdaL1; 
        else
           H(j,j+3)=lambdaL2; 
        end
    end

    z=[dd_phase1(1:nDD,i); dd_phase2(1:nDD,i)]; 
       
    x=eye(nParm)*x;
    Px=eye(nParm)*Px*eye(nParm)'+Q;
    K=Px*H'*inv(H*Px*H'+R);
    y=H*x;
    x=x+K*(z-y);
    Px=(eye(nParm)-K*H)*Px;
    
    % apply the LAMBA method in parallel
    [fix,sqr,Ps,Qz,Z,nFix,mu]=LAMBDA(x(nNonBiasParm+1:nParm),Px(nNonBiasParm+1:nParm,nNonBiasParm+1:nParm),6,'ncands',2);
    
     Fixed=fix(:,1);
     if(nFix>0)
        x_fixed(1:nNonBiasParm,i)=x(1:nNonBiasParm)-Px(1:nNonBiasParm,nNonBiasParm+1:nParm)*pinv(Px(nNonBiasParm+1:nParm,nNonBiasParm+1:nParm))*(x(nNonBiasParm+1:nParm)-fix(:,1));
        x_fixed(nNonBiasParm+1:nParm,i)=Fixed;
     else
        notFixed=notFixed+1;
        x_fixed(1:nNonBiasParm,i)=x(1:nNonBiasParm);
        x_fixed(nNonBiasParm+1:nParm,i)=Fixed;
     end
    
     nParm_=nParm;
     Px_=Px;
     x_=x;
     x_save(1:length(x_),i)=x_;
     clear Px;
     clear x;
     clear Y_sig;
     
     reference_satellite_ = reference_satellite;
end

x=x_save;
% compute absolute position estimate, float and fixed and plot comparisons
for i=1:epochs
    
    xyzEst=x(1:3,i)+tru_usrxyzO';
    estENU(:,i)=xyz2enu(xyzEst,tru_usrxyzO);
    
    rms_est(i) = sqrt(estENU(:,i)'*estENU(:,i));
    
    xyzFix=x_fixed(1:3,i)+tru_usrxyzO';
    fixENU(:,i)=xyz2enu(xyzFix,tru_usrxyzO);
    
    rms_fix(i) = sqrt(fixENU(:,i)'*fixENU(:,i));
end

sqrt(x(1:3,epochs)'*x(1:3,epochs))
figure;
plot(rms_est-enu_baseline(1));
figure;
plot(rms_fix-enu_baseline(1));




