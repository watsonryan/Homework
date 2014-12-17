clear all 
close all
clc


%% Pre-process Data

load('realDataMat')

% Define Variables 
epochs=150;%length(satTOWSEC);
lambdaL1=299792458/1575.42e6;
lambdaL2=299792458/1227.6e6;
xyz_wvmo=llh2xyz([39.64589002*pi/180,-79.96992114*pi/180,326.1251]);
xyz_pagm=llh2xyz([39.89893047*pi/180,-80.16011165*pi/180,262.7526]);
xyz_baseline=xyz_pagm;
enu_baseline = xyz2enu(xyz_baseline,xyz_wvmo)';
Sigma=3e3;
rval=.00555;
qval=0.0025;
nNonBiasParm=3;
notFixed = 0;

% Memory Allocation 

% Data Section
Sats_Used = zeros(12,288);
Used_Sats_XYZ = zeros(22,3, 288);

sd_rangeL1 = zeros(20,epochs);
sd_phaseL1 = zeros(20,epochs);
sd_rangeL2 = zeros(20,epochs);
sd_phaseL2 = zeros(20,epochs);
ref_sat_vect = zeros(1,epochs);
sd_prns = zeros(20,epochs);
dd_range1=zeros(10,epochs);
dd_range2=zeros(10,epochs);
dd_phase1=zeros(10,epochs);
dd_phase2=zeros(10,epochs);
dd_prns = zeros(20,epochs);
n_DD = zeros(1,epochs);

% KF Section 
x_fixed=zeros(20,epochs);
x_save=zeros(200,epochs);

% Data
% Range
R1_wv = wvmo.C1;
R2_wv = wvmo.P2;
R1_pa = pagw.C1;
R2_pa = pagw.P2;

% Phase
P1_wv = wvmo.PHASE1*299792458/1575.42e6;
P2_wv = wvmo.PHASE2*299792458/1227.6e6;
P1_pa = pagw.PHASE1*299792458/1575.42e6;
P2_pa = pagw.PHASE2*299792458/1227.6e6;

% Data Main Loop 
for i = 1:epochs
    
    [prn_used,nSD] = Intersection(R1_wv(:,i),R2_wv(:,i),R1_pa(:,i), R2_pa(:,i),...
                                  P1_wv(:,i),P2_wv(:,i),P1_pa(:,i),P2_pa(:,i),i);
    
    for j=1:nSD
       
        sd_prn=prn_used(j);
        
        % Calculate Singel Differences
        sd_rangeL1(j,i) = R1_pa(sd_prn,i)-R1_wv(sd_prn,i);
        sd_phaseL1(j,i) = P1_pa(sd_prn,i)-P1_wv(sd_prn,i);
        sd_rangeL2(j,i) = R2_pa(sd_prn,i)-R2_wv(sd_prn,i);
        sd_phaseL2(j,i) = P2_pa(sd_prn,i)-P2_wv(sd_prn,i);
        
        %Calculate Elevation Angle 
        sat_enu=xyz2enu([satsXYZ(i,sd_prn,1),satsXYZ(i,sd_prn,2),satsXYZ(i,sd_prn,3)], xyz_wvmo);
        el_R(j)=atan2(sat_enu(3),sqrt(sat_enu(1)^2+sat_enu(2)^2))*(180/pi);
        
    end
    
    %Matrix with Used PRN's 
    sd_prns(1:nSD,i) = prn_used;
    
    % Calculate Reference Satellites 
    
    ref_sat_candidate_vec=find(el_R==max(el_R));
    ref_sat_candidate_prn = sd_prns(ref_sat_candidate_vec,i);

    ref_sat_xyz(1,1:3,1) = [satsXYZ(i,ref_sat_candidate_prn,1)...
                           satsXYZ(i,ref_sat_candidate_prn,2)...
                           satsXYZ(i,ref_sat_candidate_prn,3)];
    ref_prn = ref_sat_candidate_prn;    
    
    %function to check elevation angle of ref satellite 
    if i == 1
        ref_prn = ref_sat_candidate_prn;
        ref_prn_ = ref_sat_candidate_prn;
    else
        [ref_prn] = ref_satellites(sd_prns(:,i-1),sd_prns(:,i),...
                                   ref_prn_,ref_prn,satsXYZ,i,xyz_wvmo);
    end
    
    nDD = nSD -1;
    reference_sat_prn(i) = ref_prn;
    
    ref_sat_vec = find(prn_used == ref_prn); 
    
    %Matrix of satellite prn numbers for double differences 
    dd_prns(1:nDD,i) = prn_used(prn_used~=ref_prn);
    
    % Calculate double differences 
    z=0;
    for j=1:nSD
        
        if(prn_used(j)~=ref_prn)
            
            dd_range1(j-z,i)=sd_rangeL1(j,i)-sd_rangeL1(ref_sat_vec,i);
            dd_range2(j-z,i)=sd_rangeL2(j,i)-sd_rangeL2(ref_sat_vec,i);
            dd_phase1(j-z,i)=sd_phaseL1(j,i)-sd_phaseL1(ref_sat_vec,i);
            dd_phase2(j-z,i)=sd_phaseL2(j,i)-sd_phaseL2(ref_sat_vec,i);
            
        else

            z=1;

         end
        
    end
    
    ref_prn_ = ref_prn;
    n_DD(i) = nDD;
    ref_sat_vect(i) = ref_sat_vec;
    clear el_R
end
%% 




%% Kalman Filter

nParm=3+2*n_DD(1);
nNonBiasParm=3;

%Initialize State Vector
x=[xyz_pagm-xyz_wvmo, dd_range1(1:n_DD(1),1)'/lambdaL1,dd_range2(1:n_DD(1),1)'/lambdaL2]';
%Initialize Error Covariance 
Px=diag([1,1,1,ones(1,n_DD(1)*2)*1000^2]);
%Initialize Noise Covariance
R=(rval^2)*(ones(n_DD(1)*2,n_DD(1)*2)+eye(n_DD(1)*2));

% Save for next iteration 
x_=x;
Px_=Px;
nParm_=nParm;

x_save(1:length(x),1) = x;
x_fixed(1:length(x),1) = x;
% Main KF Loop
for i=2:epochs
    
    nDD = n_DD(i);
    nParm=3+2*nDD;
    ref_sat_prn =reference_sat_prn(i);
   
     if i >2 && ref_sat_prn_ ~= ref_sat_prn
         
        Amb_L1 = nonzeros(dd_range1(:,i))./lambdaL1;
        Amb_L2 = nonzeros(dd_range2(:,i))./lambdaL2;
        x = [x_(1:3); Amb_L1; Amb_L2];
        
        Px = Sigma*eye(2*nDD+3);
        Px(1:3,1:3) = Px_(1:3,1:3);
        
     else
        
        [x,Px]=addRemoveBiasParm(x_(1:nParm_),Px_,nNonBiasParm,nonzeros(dd_prns(:,i-1)),nonzeros(dd_prns(:,i)), Sigma^2);
         
     end
    
     
    R=(rval^2)*(ones(nDD*2,nDD*2)+eye(nDD*2));
    Q=zeros(nParm);
    
    for j=1:3
        
        Q(j,j)=qval;
        
    end
    
    Rsat=[satsXYZ(i,ref_sat_prn,1),...
          satsXYZ(i,ref_sat_prn,2),...
          satsXYZ(i,ref_sat_prn,3)];
      
    vec2RSat=((Rsat-xyz_wvmo)/norm(Rsat-xyz_wvmo));
    uSats=[];
    
    for j=1:nDD
        
        xyzSat=[satsXYZ(i,dd_prns(j,i),1),...
                satsXYZ(i,dd_prns(j,i),2),...
                satsXYZ(i,dd_prns(j,i),3)];
            
        vecXYZSat=((xyzSat-xyz_wvmo)/norm(xyzSat-xyz_wvmo));
        uSats(j,1:3)=(vecXYZSat-vec2RSat);
        
    end
    
    uSats=[uSats;uSats]; 
        
    H=zeros(2*nDD,nParm);
    H(:,1:3)=uSats(:,1:3);
    for j=1:2*nDD
        
        if(j<=nDD)
            
            H(j,j+3)=lambdaL1;
            
        else
            
            H(j,j+3)=lambdaL2; 
            
        end
        
    end
    
    z=[dd_phase1(1:nDD,i); dd_phase2(1:nDD,i)]; 
      
    %Kalman Filter Equations
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
    
     %For next iteration 
     nParm_=nParm;
     Px_=Px;
     x_=x;
     ref_sat_prn_ = ref_sat_prn;
     ref_sat_vec_ = ref_sat_vec;
     x_save(1:length(x_),i)=x_;
     clear Px x Y_sig;

end

x=x_save;


% compute absolute position estimate, float and fixed and plot comparisons
for i=1:epochs
    
    xyzEst=x(1:3,i)+xyz_wvmo';
    estENU(:,i)=xyz2enu(xyzEst,xyz_wvmo);
    
    rms_float(i) = sqrt(estENU(:,i)'*estENU(:,i))-32481.1885;
    
    xyzFix=x_fixed(1:3,i)+xyz_wvmo';
    fixENU(:,i)=xyz2enu(xyzFix,xyz_wvmo);
    
    rms_fix(i) = sqrt(fixENU(:,i)'*fixENU(:,i))-32481.1885;
    
end

figure;
plot(rms_float);
figure;
plot(rms_fix);
    
