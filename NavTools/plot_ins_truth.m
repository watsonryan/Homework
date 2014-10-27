%
%   plot_ins_truth.m
%
clear all

load prof_1_ins_truth    % load output from GEN_INS_TRUTH.M
load flt_prof_1_dat     % load output from FLIGHT_PROFILE_GEN_1.M

close all
subplot(211)
plot(time/60,tru_lat*180/pi)
title('Baseline without sensor or initialization errors')
ylabel('latitude in degrees')
xlabel('time in minutes')

subplot(212)
plot(time/60,tru_lon*180/pi)
ylabel('longitude in degrees')
xlabel('time in minutes')

east_vel_err = tru_vel_L(:,1) - profile(:,4);
north_vel_err = tru_vel_L(:,2) - profile(:,5);
up_vel_err = tru_vel_L(:,3) - profile(:,6);

figure
subplot(311)
plot(time/60,east_vel_err,'b-')
title('Velocity Difference in Meters/Second')
ylabel('east')

subplot(312)
plot(time/60,north_vel_err,'b-')
ylabel('north')

subplot(313)
plot(time/60,up_vel_err,'b-')
ylabel('up')
xlabel('time in minutes')
