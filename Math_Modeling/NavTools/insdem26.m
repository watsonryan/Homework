%
%   insdem26.m
%
%   Simplified F-16 model
%
%   Long distance trajectory generated
%   in local-level-tangent plane (ENU)
%   coordinates.
%
%    03/03/05
%
%
clear all
close all

initpos = [0 0 10000];   % initial position
initvel = [141.4214 141.4214 0];   % initial velocity  % MUST BE CONSISTENT WITH EULER ANGLES ! ! !

phi=0*pi/180;       % initial roll angle is zero
theta=0*pi/180;     % initial pitch angle is zero
psi=45*pi/180;       % initial yaw angle

initdcm=eulr2dcm([phi theta psi]);   % Compute initial direction cosine matrix
%                                    % for aircraft attitude (nav-to-body frame)

% seg number; duration; v_final; climb_ang; turn-rate; heading_final;
% turn-delta; Tinc=0.01s
segparf16 = [5 3600  200 NaN NaN 45 NaN 1;  % Level flight for 3600 sec
             8 NaN NaN NaN 3 NaN NaN 0.01;  % Pitch and roll into a turn
             9 NaN NaN NaN 3 135 -1.913 0.005;  % Turn to heading of 135 deg
             8 NaN NaN NaN 0 NaN NaN 0.01;  % Pitch and roll back to level
             5 1  200 NaN NaN 135 NaN 0.01;  % Level flight for 1 sec
             5 3000  200 NaN NaN 135 NaN 1];  % Level flight for 3000 sec

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
title('INSDEM26:  Position Plot')
xlabel('east (meters)')
ylabel('north (meters)')
zlabel('up (meters)')
grid

figure
subplot(311)
plot(time,profile(:,4))
title('INSDEM26:  Velocity Components')
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
title('INSDEM26:  Euler Angles')
ylabel('roll angle in deg')
subplot(312)
plot(time,pitch)
ylabel('pitch angle in deg')
subplot(313)
plot(time,yaw)
ylabel('yaw angle in deg')
xlabel('run time in seconds')

