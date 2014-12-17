%  5-State GPS Unscented Kalman Filter 
%
%  The state vector is:
%   x1 = delta x position
%   x2 = delta y position
%   x3 = delta z position
%   x4 = delta clock offset
%   x5 = delta clock drift

% Authors: Ryan M. Watson,

clear all
%close all
load dataSet5.mat % static data
nom1=ones(1,length(satsXYZ));
nom1=truthXYZ(1,1).*nom1;
nom2=ones(1,length(satsXYZ));
nom2=truthXYZ(2,1).*nom2;
nom3=ones(1,length(satsXYZ));
nom3=truthXYZ(3,1).*nom3;

nomXYZ=[nom1;nom2;nom3];

prData=(2.546*prDataP1)-(1.546*prDataP2);
Ntot=length(prData);
%UKF Parameters and Weights
L = 10; % Size of state vector
alpha = 1; % Primary scaling parameter
beta = 2; % Secondary scaling parameter (Gaussian assumption)
kappa = 0; % Tertiary scaling parameter
lambda = alpha^2*(L+kappa) - L;
wm = ones(2*L + 1,1)*1/(2*(L+lambda));
wc = wm;
wm(1) = lambda/(lambda+L);
wc(1) = lambda/(lambda+L) + 1 - alpha^2 + beta;

speedoflight = 299792458; % Speed of Light (m/s)

Q=[.01; .01; .01; .01; .01];
Q = diag(Q);
sigmaz = 2; % standard deviation of pseudorange measurement
% noise in meters

x_pre_ukf = zeros(10,Ntot);  % initial state prediction
P_pre = 50*eye(5);   % initial prediction error covariance matrix
P=P_pre;
x_ukf=x_pre_ukf;

Pvec=[];
prevsm=[];
prevadr=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%   Non-Additive Noise UKF   %%%%%%%%%%%%%%%%%%%%%%%
    
Ntot=length(prData);
for i = 1:Ntot
    if i==1
        estusr = olspos(prData(1:nSat(i),i),satsXYZ(1:nSat(i),:,i));
%         disp(estusr)
        x_ukf(1,1)=estusr(1);
        x_ukf(2,1)=estusr(2);
        x_ukf(3,1)=estusr(3);
        x_ukf(4,1)=estusr(4);
    else
        
        R = sigmaz*eye(nSat(i));   % form the measurement error covariance matrix
        
        [prsm,Pvec,prevadr]=compkalm(prData(1:nSat(i),i),prRateL1(1:nSat(i),i),svID(1:nSat(i),i),1,.01,2,Pvec,prevsm,prevadr)
        Smooth(:,i)=prsm;
        Pvec=Pvec;
        prsm=prsm;
        prevadr=prevadr;
        
        [sP,dummy] = chol(P,'lower'); % Calculate square root of error covariance
        sQ = sqrtm(Q);
        sPQ = blkdiag(sP,sQ);
        % chi_p = "chi previous" = chi(k-1)
        chi_p = [x_ukf(:,i-1), x_ukf(:,i-1)*ones(1,L)+sqrt(L+lambda)*sPQ,x_ukf(:,i-1)*ones(1,L)-sqrt(L+lambda)*sPQ];
        
        
        chi_p_Smooth = [Smooth(:,i-1), Smooth(:,i-1)*ones(1,L)+sqrt(L+lambda)*sPQ,Smooth(:,i-1)*ones(1,L)-sqrt(L+lambda)*sPQ];
        
        for j = 1:2*L+1
            %State Update
            % position random walk
            chi_p(1:3,j) = chi_p(1:3,j) + chi_p(6:8,j);
            % white noise drift
            chi_p(5,j) = chi_p(5,j)+chi_p(10,j);
            % bias = bias+drift
            chi_p(4,j) = chi_p(4,j)+chi_p(9,j) + chi_p(5,j);
        end
        
        x_ukf(1:5,i) = chi_p(1:5,:)*wm; % Calculate mean of predicted state
        % Calculate covariance of predicted stats
        P_m = zeros(5,5);
        for k = 1:2*L+1
            P_m = P_m + wc(k)*(chi_p(1:5,k) - x_ukf(1:5,i))*(chi_p(1:5,k) - x_ukf(1:5,i))';
        end
        
        
        for j = 1:2*L+1
            %State Update
            % position random walk
            chi_p_Smooth(1:3,j) = chi_p_Smooth(1:3,j) + chi_p_Smooth(6:8,j);
            % white noise drift
            chi_p_Smooth(5,j) = chi_p_Smooth(5,j)+chi_p_Smooth(10,j);
            % bias = bias+drift
            chi_p_Smooth(4,j) = chi_p_Smooth(4,j)+chi_p_Smooth(9,j) + chi_p_Smooth(5,j);
        end
        
        Smooth(1:5,i) = chi_p_Smooth(1:5,:)*wm; % Calculate mean of predicted state
        % Calculate covariance of predicted stats
        P_m_Smooth = zeros(5,5);
        for k = 1:2*L+1
            P_m_Smooth = P_m_Smooth + wc(k)*(chi_p_Smooth(1:5,k) - Smooth(1:5,i))*(chi_p_Smooth(1:5,k) - Smooth(1:5,i))';
        end
        
        L_1 = 5+nSat(i); % Size of state vector
        lambda_1 = alpha^2*(L_1+kappa) - L_1;
        wm_1 = ones(2*L_1 + 1,1)*1/(2*(L_1+lambda_1));
        wc_1 = wm_1;
        wm_1(1) = lambda_1/(lambda_1+L_1);
        wc_1(1) = lambda_1/(lambda_1+L_1) + 1 - alpha^2 + beta;
        
        
        [sP_m, dum] = chol(P_m,'lower'); % Calculate square root of error covariance
        sR = sqrtm(R);
        sPR = blkdiag(sP_m,sR);
       
        x_m = zeros(L_1,1);
        x_m(1:5,1)=x_ukf(1:5,i);
        chi_m = [x_m, x_m*ones(1,L_1)+sqrt(L_1+lambda_1)*sPR,x_m*ones(1,L_1)-sqrt(L_1+lambda_1)*sPR];
        
        x_m_Smooth = zeros(L_1,1);
        x_m_Smooth(1:5,1)=Smooth(1:5,i);
        chi_m = [x_m_Smooth, x_m_Smooth*ones(1,L_1)+sqrt(L_1+lambda_1)*sPR,x_m_Smooth*ones(1,L_1)-sqrt(L_1+lambda_1)*sPR];
        
            for j=1:nSat(i)
                TropoDel(j,i)= tropocorr(satsXYZ(j,:,i),x_ukf(1:3,i));
                TropoDelSmooth(j,i)=tropocorr(satsXYZ(j,:,i),Smooth(:,i));
            end
 
    
        Y_sig = zeros(nSat(i),1+2*L_1);
        for k = 1:2*L_1+1
            for l=1:nSat(i)
                Y_sig(l,k) = norm(satsXYZ(l,:,i) - chi_m(1:3,k)')+chi_m(4,k)+chi_m(5+l,k)+TropoDel(l,i);
            end
        end
        
        
         Y_sig_Smooth = zeros(nSat(i),1+2*L_1);
        for k = 1:2*L_1+1
            for l=1:nSat(i)
                Y_sig_Smooth(l,k) = norm(satsXYZ(l,:,i) - Smooth(:,k)')+chi_m(4,k)+chi_m(5+l,k)+TropoDelSmooth(l,i);
            end
        end
        
        y_m(1:nSat(i),1) = Y_sig(1:nSat(i),:)*wm_1;
        Pyy = zeros(nSat(i));
        for k = 1:2*L_1+1
            Pyy = Pyy + wc_1(k)*(Y_sig(1:nSat(i),k) - y_m(1:nSat(i),1))*(Y_sig(1:nSat(i),k) - y_m(1:nSat(i),1))';
        end
        
         y_m_Smooth(1:nSat(i),1) = Y_sig_Smooth(1:nSat(i),:)*wm_1;
        Pyy_Smooth = zeros(nSat(i));
        for k = 1:2*L_1+1
            Pyy_Smooth = Pyy + wc_1(k)*(Y_sig_Smooth(1:nSat(i),k) - y_m_Smooth(1:nSat(i),1))*(Y_sig_Smooth(1:nSat(i),k) - y_m_Smooth(1:nSat(i),1))';
        end
        
        Pxy = zeros(5,nSat(i));
        for k = 1:2*L_1+1
            Pxy = Pxy + wc_1(k)*(chi_m(1:5,k)-x_m(1:5,1))*(Y_sig(1:nSat(i),k) - y_m(1:nSat(i),1))';
        end
        
        Pxy_Smooth = zeros(5,nSat(i));
        for k = 1:2*L_1+1
            Pxy_Smooth = Pxy_Smooth + wc_1(k)*(chi_m(1:5,k)-x_m(1:5,1))*(Y_sig_Smooth(1:nSat(i),k) - y_m_Smooth(1:nSat(i),1))';
        end
        
        K = Pxy/Pyy;
        prz=prData(1:nSat(i),i);
        x_ukf(1:5,i) = x_m(1:5,1)+K*(prz(1:nSat(i))-y_m(1:nSat(i)));
        P = P_m - K*Pyy*K';
        
        K_Smooth = Pxy_Smooth/Pyy_Smooth;
        prz_Smooth=Smooth(1:nSat(i),i);
        x_ukf_Smooth(1:5,i) = x_m(1:5,1)+K_Smooth*(prz_Smooth(1:nSat(i))-y_m_Smooth(1:nSat(i)));
        P_Smooth = P_m - K_Smooth*Pyy_Smooth*K_Smooth';
        
    end
end



for i = 1:length(satsXYZ)
ENU_UKF = xyz2enu(x_ukf(1:3,i),nomXYZ(:,i));
ENU_UKF_Smooth = xyz2enu(Smooth(:,i),nomXYZ(:,i));
ENU_Truth = xyz2enu(truthXYZ(:,i),nomXYZ(:,i));
ERROR(:,i)=ENU_UKF-ENU_Truth;
ERROR_SMOOTH(:,i)=ENU_UKF_Smooth-ENU_Truth;
Clock = x_ukf(4,:)-truthClockBias;
end

Clock(1,1)=0; ERROR(3,1)=0;

Mean_Error_East=mean(ERROR(1,:))
Mean_Error_North=mean(ERROR(2,:))
Mean_Error_UP=mean(ERROR(3,:))
Mean_Clock_Error=mean(Clock)

subplot(4,1,1)
plot(time,ERROR(1,:),'LineWidth',2, 'Color', [0,0,0])
hold on
plot(time,ERROR_SMOOTH(1,:),'LineWidth',2, 'Color', [0,0,0])
title('East')
axis tight
subplot(4,1,2)
plot(time,ERROR(2,:),'LineWidth',2, 'Color', [0,0,0])
hold on
plot(time,ERROR_SMOOTH(2,:),'LineWidth',2, 'Color', [0,0,0])
title('North')
axis tight
subplot(4,1,3)
plot(time,ERROR(3,:),'LineWidth',2, 'Color', [0,0,0])
hold on
plot(time,ERROR_SMOOTH(3,:),'LineWidth',2, 'Color', [0,0,0])
title('UP')
axis tight
subplot(4,1,4)
plot(time,Clock,'LineWidth',2, 'Color', [0,0,0])
title('Clock Bias')
axis tight



