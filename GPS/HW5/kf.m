clear all; close all; clc;
load realDataMat.mat
	
xyz_wvmo=[39.64589002*pi/180,-79.96992114*pi/180,326.1251];%llh2xyz([39.64589002*pi/180,-79.96992114*pi/180,326.1251]);
idx=findHighestElevationSat(satsXYZ, xyz_wvmo);
[dd_range1, dd_range2, dd_phase1, dd_phase2]=findDD(wvmo,pagw,idx);

%% UKF Parameters and Weights
L = 3+2*numDD+3; % Size of augmented state vector pos + ambiguities + process noise 
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

%% UKF Initializations
numSats=-------------------------------
numEpochs=--------------------------------

x=zeros(3+2*numSats,1000);
x(:,1)=[0,0,0, dd_range1(1,:)/299792458/1575.42e6,ddrange2(1,:)/299792458/1227.6e6]';
x_fixed=zeros(2*numSats,numEpochs);

notFixed=0;

Px=diag([1,1,1,ones(1,numDD*2)*15^2]);
Q=eye(3)*2;
R=(.02^2)*(ones(numDD*2,numDD*2)+eye(numDD*2));

[sR,dum]=sqrtm(R);  % only needed once, so dont do in loop

L_meas = 3+4*numDD; % size of state vector
lambda_meas = alpha^2*(L_meas+kappa) - L_meas;
wm_meas = ones(2*L_meas + 1,1)*1/(2*(L_meas+lambda_meas));
wc_meas = wm_meas;
wm_meas(1) = lambda_meas/(lambda_meas+L_meas);
wc_meas(1) = lambda_meas/(lambda_meas+L_meas) + 1 - alpha^2 + beta;

for i=2:numEpochs
    % create the sigma-points
    % copy the state vector L times
    XA = zeros(L,L);
    for j=1:L;
        XA(1:nParm,j)=x(1:nParm,i-1);
    end
    Pk=blkdiag(Px,Q);
    [rootMat ,dum]=sqrtm(Pk);
    sP=sqrt(L+lambda)*rootMat;
    
    % add and subtract sqrt(P) to produce the sigma-points
    x_nom0=[x(1:nParm,i-1); zeros(length(Q),1)];
    X_sig=[x_nom0 XA+sP XA-sP];
    
    % run the state prediction with each set of sigma-points
    for j=1:2*L+1;
        % our process model is the identity plus uncertainty
        X_sig(1:3,j)=X_sig(1:3,j)+X_sig(nParm+1:nParm+3,j);
        
    end
    
    % calculate the predicted state vector as a weighted average of
    % the predicted sigma points
    x(1:nParm,i)=X_sig(1:nParm,:)*wm;
    
    % computed the predicted error-covariance using the transformed sigma
    % points
    Px=zeros(nParm,nParm);
    for j=1:2*L+1;
        Px=Px+wc(j)*((X_sig(1:nParm,j)-x(:,i))*(X_sig(1:nParm,j)-x(:,i))');
    end
    
    % determine the unit vector to the reference satellite ( only needed
    % once per epoch )
    Rsat=[sats(find(id(i,:)==refSat),1,i) sats(find(id(i,:)==refSat),2,i) sats(find(id(i,:)==refSat),3,i)];
    vec2RSat=((Rsat-tru_usrxyzO)/norm(Rsat-tru_usrxyzO));
    
    % measurement vector only carrier phase double differences
    z=[ddPhaseL1(i,1:nSat(i)-1) ddPhaseL2(i,1:nSat(i)-1) R_uwb(i)]';
    
    
    % show a sequential update Kalman Filter ( useful when more or less
    % sats included )
    lambdaL1=299792458/1575.42e6;
    lambdaL2=299792458/1227.6e6;
    for ii=1:nObs;
        
        
        % create the measurement update sigma points using the
        % predicted state and error covariance
        [sPx,dum]=sqrtm(real(Px));
        Pmk=blkdiag(sPx,sR);
        
        sP_meas = sqrt(L_meas+lambda_meas)*Pmk;
        XA_meas = zeros(L_meas,L_meas);
        x_nom0=[x(:,i); zeros(nObs,1)];
        % copy the predicted state vector
        for j=1:L_meas;
            XA_meas(:,j)=x_nom0;
        end
        % add and substract the sqrt(P)
        X_sig_meas=[x_nom0 XA_meas+sP_meas XA_meas-sP_meas];
        
        % use the predicted baseline and ambiguities to compute
        % what expected phase measurement is
        for j=1:2*L_meas+1;
            
            Ind=mod(ii-1,nSat(i)-1)+1;
            Ssat=[sats(find(id(i,:)==idSVdd(i,Ind)),1,i) sats(find(id(i,:)==idSVdd(i,Ind)),2,i) sats(find(id(i,:)==idSVdd(i,Ind)),3,i)];
            vec2Sat=((Ssat-tru_usrxyzO)/norm(Ssat-tru_usrxyzO));
            
            % first set of measurements are L!
            if ii<=numDD
                % L1 computed Phase equals -difference unit vectors * predicted
                % baseline + predicted ambiguity* lamndaL1
                Y_sig(ii,j)=-(vec2Sat-vec2RSat)*X_sig_meas(1:3,j)+X_sig_meas(nParm+ii,j)+X_sig_meas(ii+nNonBiasParm,j)*lambdaL1;
            else
                % L2 computed Phase equals -difference unit vectors * predicted
                % baseline + predicted ambiguity* lamndaL2
                Y_sig(ii,j)=-(vec2Sat-vec2RSat)*X_sig_meas(1:3,j)+X_sig_meas(nParm+ii,j)+X_sig_meas(ii+nNonBiasParm,j)*lambdaL2;
            end
            
        end
        
        % because we are doing sequential updates, we are now dealing with
        % scalars and vectors
        y=Y_sig(ii,:)*wm_meas;
        Py=0;
        for j=1:2*L_meas+1;
            Py=Py+wc_meas(j)*((Y_sig(ii,j)-y)*(Y_sig(ii,j)-y)');
        end
        Pxy=zeros(nParm,1);
        for j=1:2*L_meas+1;
            Pxy=Pxy+wc_meas(j)*((X_sig_meas(1:nParm,j)-x(1:nParm,i))*(Y_sig(ii,j)-y)');
        end
        K=Pxy/Py;
        innovationResidual(i,ii)=z(ii)-y;
        
        Px_=Px;
        
        x(1:nParm,i)=x(1:nParm,i)+K*(innovationResidual(i,ii));
        
        Px=Px-K*Py*K';
    end % end sequential measurement update loop
    
    % apply the LAMBA method in parallel
    [fix,sqr,Ps,Qz,Z,nFix,mu]=LAMBDA(x(nNonBiasParm+1:nParm,i),Px(nNonBiasParm+1:nParm,nNonBiasParm+1:nParm),6,'ncands',2);
    
    Fixed(:,i)=fix(:,1);
    if(nFix>0)
    x_fixed(1:nNonBiasParm,i)=x(1:nNonBiasParm,i)-Px(1:nNonBiasParm,nNonBiasParm+1:nParm)*pinv(Px(nNonBiasParm+1:nParm,nNonBiasParm+1:nParm))*(x(nNonBiasParm+1:nParm,i)-fix(:,1));
    x_fixed(nNonBiasParm+1:nParm,i)=Fixed(:,i);
    else
    notFixed=notFixed+1;
    x_fixed(1:nNonBiasParm,i)=x(1:nNonBiasParm,i);
    x_fixed(nNonBiasParm+1:nParm,i)=Fixed(:,i);
    end
   
end

% compute absolute position estimate, float and fixed and plot comparisons
for i=1:numEpochs
    
    xyzEst=deltaXYZR+x(1:3,i);
    estENU(:,i)=xyz2enu(xyzEst,tru_usrxyzO);
    
    xyzFix=deltaXYZR+x_fixed(1:3,i);
    fixENU(:,i)=xyz2enu(xyzFix,tru_usrxyzO);
    
end


figure,
plot(trueENU(1,:),trueENU(2,:),'lineWidth',3)
hold on
plot(estENU(1,:),estENU(2,:),'r--','lineWidth',3)
plot(fixENU(1,:),fixENU(2,:),'g','lineWidth',3)
xlabel('East (m)')
ylabel('North (m)')
legend('Float','Fixed')
grid on
grid on


figure,plot((estENU-trueENU)')
hold on
grid on
plot((fixENU-trueENU)','--','lineWidth',3)
