%% Andrew Saiko - MAE 593 - GPS - Homework #2
format long g
clear
clc
LOAD=menu('Select which data set to load for LLS Model','Data Set 3','Data Set 4');
if LOAD==1
load('f:\MAE 593\HW2\dataSet3.mat')
end
if LOAD==2
load('f:\MAE 593\HW2\dataSet4.mat')
end
%% Creates a 14x3 matrix with the nominal position
nom_XYZ=ones(14,3);
nom_XYZ(:,1)=nomXYZ(1);
nom_XYZ(:,2)=nomXYZ(2);
nom_XYZ(:,3)=nomXYZ(3);
nomCLOCK=clockBiasNom;
%% Establishes Constants needed for Trop Delay
P0=1013; %Total pressure [mbar]
T0=288.15; %Temperature [K]
e0=12.8; %Partial pressure due to water [mbar]
%% Begins iterations
for i=1:2400;
    %% Truth Inputs for each iteration, i
    ClockBias_Truth = truthClockBias(i);
    XYZ_Truth = truthXYZ(:,i);
    %% Transforming XYZ -> ENU
    [enu] = xyz2enu(XYZ_Truth,nomXYZ);
    enu_TRUTH(:,i) = enu;
    %% Transforming XYZ -> Lat, Lon, Height
    [llh] = xyz2llh(XYZ_Truth);
    llh_TRUTH(i,:) = llh;
    %% Giving Inputs and Calculating XYZ/CLOCK Estimates
    numSat=nSat(i);
    pr_OBS = prDataP1(1:nSat(i),i);
    sats_XYZ = satsXYZ(1:nSat(i),1:3,i);
    c = 299792458; % m/s
    [Clock_Bias_Estimate,XYZ_Estimate] = GPS(nomXYZ,nomCLOCK,pr_OBS,sats_XYZ,numSat);
    %% Storing Function Estimates in a matrix
    ClockBias_Est(i,1) = Clock_Bias_Estimate;
    XYZ_Est(i,:) = XYZ_Estimate;
    %% Calculating Errors and storing in a matrix
    I_Error_XYZ(i,1) = norm(nomXYZ'-truthXYZ(:,i));
    F_Error_XYZ(i,1) = norm(XYZ_Est(i,:)'-truthXYZ(:,i));
    I_Error_CLOCK(i,1) = norm(nomCLOCK-truthClockBias(1,i));
    F_Error_CLOCK(i,1) = norm(ClockBias_Est(i,1)-truthClockBias(1,i));
    %% Transforming Estimated XYZ -> ENU
    [enu] = xyz2enu(XYZ_Estimate,nomXYZ);
    enu_Est(:,i) = enu;
    %% Transforming Estimated XYZ -> Lat, Lon, Height
    [llh] = xyz2llh(XYZ_Estimate);
    llh_Est(:,i) = enu;
    %% Calculated ENU Errors
    I_Error_XYZ_ENU(i,1) = norm(nomXYZ'-enu_TRUTH(:,i));
    F_Error_XYZ_ENU(i,1) = norm(enu_Est(:,i)-enu_TRUTH(:,i));
end
%% Calculates Error in ENU
E_east=enu_Est(1,:)-enu_TRUTH(1,:);
E_north=enu_Est(2,:)-enu_TRUTH(2,:);
E_up=enu_Est(3,:)-enu_TRUTH(3,:);
E_clock=ClockBias_Est-truthClockBias';
AVGe=mean(E_east);
AVGn=mean(E_north);
AVGu=mean(E_up);
AVGc=mean(E_clock);
AVG_E=[AVGe AVGn AVGu];
%% Calculating RMS Errors
RMSe=rms(E_east);
RMSn=rms(E_north);
RMSu=rms(E_up);
RMS_ENU=[RMSe RMSn RMSu];
RMS_3D=rms(F_Error_XYZ_ENU);
RMS_clock=rms(E_clock);
%% Plot the East, North, Up Error in Subplots
Input=menu('Select the Graph to Display','Non-Tropospheric ENU Model','Comparison of ENU Models','Non-Tropospheric Clock Model','Comparison of Clock Models');

if Input==1
    subplot(3,1,1)
    plot(time,E_east)
    xlabel('Time [s]');
    ylabel('Error in East direction [m]');
    legend('No Tropospheric Model');
    
    subplot(3,1,2)
    plot(time,E_north)
    xlabel('Time [s]');
    ylabel('Error in North direction [m]');
    legend('No Tropospheric Model');
    
    subplot(3,1,3)
    plot(time,E_up)
    xlabel('Time [s]');
    ylabel('Error in Up direction [m]');
    legend('No Tropospheric Model');
    hold on
end

if Input==2
    hold on
    subplot(3,1,1);
    plot(time,E_east);
    xlabel('Time [s]');
    ylabel('Error in East direction [m]');
    legend('No Trop Model');
    
    hold on
    subplot(3,1,2);
    plot(time,E_north);
    xlabel('Time [s]');
    ylabel('Error in North direction [m]');
    legend('No Trop Model');
    
    hold on
    subplot(3,1,3);
    plot(time,E_up);
    xlabel('Time [s]');
    ylabel('Error in Up direction [m]');
    legend('No Trop Model');
    hold on
    %% Starts Tropospheric Script !!!
    Run_GPS_Trop
end
%% Plot the 3D Norm error results - ENU
% plot(time,F_Error_XYZ_ENU,time,F_Error_CLOCK)
% legend('Position Error','Clock Error')
%% Plot of the 3D Position Error
%hold on
% plot(time,F_Error_XYZ_ENU)
% legend('Position Error')
%% Plot of the Clock Error
if Input==3
    plot(time,E_clock)
    legend('Clock Error w/o Trop Model')
    xlabel('Time [s]')
    ylabel('Clock Error [m]')
end
if Input==4
    hold on
    plot(time,E_clock)
    legend('Clock Error w/o Trop Model')
    xlabel('Time [s]')
    ylabel('Clock Error [m]')
    %% Starts Tropospheric Script !!!
    Run_GPS_Trop
end
%% Calculating Difference Between Estimate and Truth ENU
% Diff_ENU = enu_TRUTH - enu_Est;
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