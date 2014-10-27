%
%  decim01.m
%
%  Program to decimate multi-time-step
%  processed vehicle trajectory data 
%  down to a uniform sampling rate of 1 Hz
%
clear all
close all

load doglegall

k = 1;
tru_time(k) = time(1);
tru_pos_L(k,:) = profile(1,1:3);
tru_vel_L(k,:) = profile(1,4:6);

for i = 2:length(time),
    if ( rem(time(i),1)<0.0005 ) | ( (time(i)-tru_time(k))>0.9995 ),
        k = k + 1;
        tru_time(k) = time(i);
        tru_pos_L(k,:) = profile(i,1:3);
        tru_vel_L(k,:) = profile(i,4:6);
    end
end

% save dogleg1hz tru_time tru_pos_L tru_vel_L -V6
save dogleg1hz tru_time tru_pos_L tru_vel_L

