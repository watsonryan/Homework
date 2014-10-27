%
%   insdem24.m
%
%   Same as INSDEM22 but a faster sampling rate during the turn
%
%   04/03/05
%
%
clear all
close all

initpos = [0 0 0];   % initial position
initvel = [0 0 0];   % initial velocity

phi=0*pi/180;       % initial roll angle is zero
theta=0*pi/180;     % initial pitch angle is zero
psi=90*pi/180;       % initial yaw angle

initdcm=eulr2dcm([phi theta psi]);   % Compute initial direction cosine matrix
%                                    % for aircraft attitude (nav-to-body frame)

% seg number; duration; v_final; climb_ang; turn-rate; heading_final;
% turn-delta; Tinc=0.01s
segparf16 = [5 20  200 NaN NaN 90 NaN 0.01;  % Level accel for 20 sec up to 200 meters/sec
             6 NaN NaN 5   NaN 90 NaN 0.01; % Pitch up transition to 5 deg climb
             7 1  NaN 5   NaN 90 NaN 0.01; % Climb for 1 sec
             7 300  NaN 5   NaN 90 NaN 1; % Climb for 300 sec
             6 NaN NaN 0   NaN 90 NaN 0.01; % Pitch down transition to level flight
             5 1  200 NaN NaN 90 NaN 0.01;  % Level flight for 1 sec
             5 60  200 NaN NaN 90 NaN 1;  % Level flight for 60 sec
             8 NaN NaN NaN -3 NaN NaN 0.005;  % Pitch and roll into a turn
             9 NaN NaN NaN -3 270 -1.59 0.005;  % Turn to heading of 270 deg
             8 NaN NaN NaN 0 NaN NaN 0.005;  % Pitch and roll back to level
             5 1  200 NaN NaN 270 NaN 0.005;  % Level flight for 1 sec
             5 60  200 NaN NaN 270 NaN 1;  % Level flight for 60 sec
             6 NaN NaN -15   NaN 270 NaN 0.01; % Pitch down transition to 15 deg descent
             7 1  NaN -15   NaN 270 NaN 0.01; % Descend for 1 sec
             7 30  NaN -15   NaN 270 NaN 1; % Descend for 30 sec
             6 NaN NaN 0   NaN 270 NaN 0.01; % Pitch up transition to level flight
             5 1  200 NaN NaN 270 NaN 0.01;  % Level flight for 1 sec
             5 90  200 NaN NaN 270 NaN 1];  % Level flight for 90 sec

profile = progenf16(initpos,initvel,initdcm,segparf16);
time = profile(:,19);   % run time in seconds

h = waitbar(0,'Calculating Euler Angles');
ntot = size(profile,1);
for k = 1:ntot,
    dcmnb=[profile(k,10:12); profile(k,13:15); profile(k,16:18)];
    dcmbn=dcmnb';
    eulv=dcm2eulr(dcmbn);
    roll(k) = eulv(1)*180/pi;
    pitch(k) = eulv(2)*180/pi;
    yaw(k) = eulv(3)*180/pi; 
    if yaw(k) < 0, yaw(k)=yaw(k)+360; end
    waitbar(k/ntot,h)
end
close(h)


figure
plot3(profile(:,1),profile(:,2),profile(:,3))
axis equal
title('INSDEM24: Flight Path Generated By PROGENF16')
xlabel('east (meters)')
ylabel('north (meters)')
zlabel('up (meters)')
grid

figure
subplot(311)
plot(time,profile(:,4))
title('INSDEM24: Velocity Components')
ylabel('east vel in m/s')
subplot(312)
plot(time,profile(:,5))
ylabel('north vel in m/s')
subplot(313)
plot(time,profile(:,6))
ylabel('vertical vel in m/s')
xlabel('run time in seconds')

figure
subplot(311)
plot(time,roll)
title('INSDEM24: Euler Angles')
ylabel('roll angle in deg')
subplot(312)
plot(time,pitch)
ylabel('pitch angle in deg')
subplot(313)
plot(time,yaw)
ylabel('yaw angle in deg')
xlabel('run time in seconds')

