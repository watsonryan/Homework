%
close all
%
figure
plot(time,profile(:,4),'-',time,profile(:,5),'-',time,profile(:,6),'-')
title('Velocity Components')
ylabel('velocity in m/s')
xlabel('run time in seconds')
text(57,-6,'east')
text(57,6,'north')
text(25,3,'up')

figure
plot(time,profile(:,7))
title('East Acceleration')
ylabel('acceleration in m/s/s')
xlabel('run time in seconds')
axis([0 80 -3 1])

figure
plot(time,profile(:,8))
title('North Acceleration')
ylabel('acceleration in m/s/s')
xlabel('run time in seconds')
axis([0 80 -3 1])

figure
plot(time,profile(:,9))
title('Up Acceleration')
ylabel('acceleration in m/s/s')
xlabel('run time in seconds')
axis([0 80 -3 1])

figure
plot(time,pitch)
title('Pitch Angle')
ylabel('pitch in degrees')
xlabel('run time in seconds')

figure
plot(time,roll)
title('Roll Angle')
ylabel('roll in degrees')
xlabel('run time in seconds')

figure
plot(time,yaw)
title('Yaw Angle')
ylabel('yaw in degrees')
xlabel('run time in seconds')
