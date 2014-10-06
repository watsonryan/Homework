% Script File for Function - GPS(nomXYZ,nomCLOCK,pr_OBS,satsXYZ,nSat)
%% Andrew Saiko - MAE 593 - GPS - Homework #1
format long g
clear
clc
load('f:\MAE 593\HW1\dataSet1.mat')
% This loads the data provided for the homework
%% Creates output variable matrix for each iteration to populate
i=1;
c=1;
x=1;
nom_XYZ=ones(14,3);
nom_XYZ(:,1)=nomXYZ(1);
nom_XYZ(:,2)=nomXYZ(2);
nom_XYZ(:,3)=nomXYZ(3);
ClockBias_Est = zeros(2400,1);
XYZ_Est = zeros(2400,3);
nomCLOCK=clockBiasNom;
I_Error_XYZ = zeros(2400,1);
F_Error_XYZ = zeros(2400,1);
I_Error_CLOCK = zeros(2400,1);
F_Error_CLOCK = zeros(2400,1);
enu_Truth = zeros(14,3,2400);
enu_Est = zeros(3,2400);
enu_TRUTH = zeros(3,2400);
PDOP = zeros(2400,1);
TDOP = zeros(2400,1);
GDOP = zeros(2400,1);
sinEL= zeros(2400,1);
tanAZ= zeros(2400,1);
WW=zeros(nSat(i));
W_Store6=zeros(6,6,456);
W_Store7=zeros(7,7,1944);

while i<2401;
    %% Transforming XYZ_Truth to ENU_Truth to establish weighted matrix
%     x=1;
%     while x<8
%         [enu, R] = xyz2enu(satsXYZ(x,:,i),nomXYZ);
%         enu_Truth(x,:,i) = enu;
%         sinEL(x,1) = enu_Truth(x,3,i)/sqrt(enu_Truth(x,1,i)^2+enu_Truth(x,2,i)^2+enu_Truth(x,3,i)^2);
%         tanAZ = enu_Truth(x,1,i)/enu_Truth(x,2,i);
%         WW(x,x) = sinEL(x,1);
%         x=x+1;
%     end
%     W=(WW(1:nSat(i),1:nSat(i)));
%     if nSat(i)==6;
%     W_Store6(:,:,i)=W;
%     elseif nSat(i)==7;
%         W_Store7(:,:,c)=W;
%         c=c+1;
%     end
    %% Truth Inputs for each iteration, i
    ClockBias_Truth = truthClockBias(i);
    XYZ_Truth = truthXYZ(:,i);
    %% Giving Inputs and Calculating XYZ/CLOCK Estimates
    numSat=nSat(i);
    pr_OBS = prData(1:nSat(i),i);
    sats_XYZ = satsXYZ(1:nSat(i),1:3,i);
    c = 299792458; % m/s
    % [XYZ_Estimate] = GPS_XYZ(nomXYZ,nomCLOCK,pr_OBS,sats_XYZ,numSat);
    % [Clock_Bias_Estimate] = GPS_CLOCK(nomXYZ,nomCLOCK,pr_OBS,sats_XYZ,numSat);
    [Clock_Bias_Estimate,XYZ_Estimate] = GPS_Updated(nomXYZ,nomCLOCK,pr_OBS,sats_XYZ,numSat,sigma_URE,W);
    %% Storing Function Estimates in a matrix
    ClockBias_Est(i,1) = Clock_Bias_Estimate;
    XYZ_Est(i,:) = XYZ_Estimate;
    PDOP(i,1) = pdop;
    TDOP(i,1) = tdop;
    GDOP(i,1) = gdop;
    %% Calculating Errors and storing in a matrix
    I_Error_XYZ(i,1) = norm(nomXYZ'-truthXYZ(:,i));
    F_Error_XYZ(i,1) = norm(XYZ_Est(i,:)'-truthXYZ(:,1));
    I_Error_CLOCK(i,1) = norm(nomCLOCK-truthClockBias(1,i));
    F_Error_CLOCK(i,1) = norm(ClockBias_Est(i,1)-truthClockBias(1,i));
    %% Transforming XYZ -> ENU
    [enu] = xyz2enu(XYZ_Estimate,nomXYZ);
    enu_Est(:,i) = enu;
    [enu] = xyz2enu(XYZ_Truth,nomXYZ);
    enu_TRUTH(:,i) = enu;
    %% Calculated ENU Errors
    I_Error_XYZ_ENU(i,1) = norm(nomXYZ'-enu_TRUTH(:,i));
    F_Error_XYZ_ENU(i,1) = norm(enu_Est(:,i)-enu_TRUTH(:,1));
    %% Iteration Counter
    i=i+1;
    %% Iterative Method
    % if i>1
    %     nomXYZ = XYZ_Estimate;
    %     nomCLOCK = Clock_Bias_Estimate;
    % end
end

%% Plot PDOP, TDOP, and GDOP
% figure
% plot (time, PDOP, time, TDOP, time, GDOP)
% xlabel('Time [s]');
% ylabel('Dilution of Precision');
% legend('Position','Time','Geometric')
%% Plot the error results - XYZ
figure
plot(time,F_Error_XYZ,time,F_Error_CLOCK)
legend('Position Error','Clock Error')
xlabel('Time [s]')
ylabel('Error [m]')
%% Plot the error results - ENU
figure
plot(time,F_Error_XYZ_ENU)
legend('Position Error')
xlabel('Time [s]')
ylabel('ENU Error [m]')
%% Calculating Difference between XYZ and ENU
figure
Diff = F_Error_XYZ - F_Error_XYZ_ENU;
plot(time,Diff)
legend('Difference in Position Error')
xlabel('Time [s]')
ylabel('Difference in Error between XYZ and ENU')
%% Calculate Difference Between Estimate and Truth
% Difference=truthXYZ'-XYZ_Est;
% Xavg=mean(Difference(:,1));
% Yavg=mean(Difference(:,2));
% Zavg=mean(Difference(:,3));
% AVG_Error=[Xavg,Yavg,Zavg]
% grid on
% scatter3(Difference(:,1),Difference(:,2),Difference(:,3),'r.')
% xlabel=('Error in X-Direction [m]');
% ylabel=('Error in Y-Direction [m]');
% zlabel=('Error in Z-Direction [m]');
% hold on
% scatter3(Xavg,Yavg,Zavg,'b*')
% legend('Error','Average')
% hold off
%% Calculating Difference Between Estimate and Truth ENU
% Diff_ENU = enu_Truth - enu_Est;
% Diff_ENU = Diff_ENU';
% scatter3(Diff_ENU(:,1),Diff_ENU(:,2),Diff_ENU(:,3),'r.')
% Xavg=mean(Diff_ENU(:,1));
% Yavg=mean(Diff_ENU(:,2));
% Zavg=mean(Diff_ENU(:,3));
% AVG_Error=[Xavg,Yavg,Zavg]
% grid on
% xlabel=('Error in X-Direction [m]');
% ylabel=('Error in Y-Direction [m]');
% zlabel=('Error in Z-Direction [m]');
% hold on
% scatter3(Xavg,Yavg,Zavg,'b*')
% legend('Error','Average')