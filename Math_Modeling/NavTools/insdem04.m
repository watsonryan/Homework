%
%   insdem04: 3-d navigation     QUATERNION UPDATING
%
%     Vehicle starts from a stationary position, undergoes
%     constant acceleration for 10 seconds, then . . .
%
%     - Non-rotating, flat earth without gravity assumed
%
%     - Locally-level, east-north-up coordinate system
%       (x: east, y: north, z: up)
%     - Heading (i.e., body attitude) is measured clockwise
%       from north
%
%     The portions of the simulation not dealing directly
%     with the quaternion are commented/documented in
%     insdem02.m

clear all
insdem01
close all

npts = max(size(time));

fprintf(1,' Creating delta-v and delta-theta profiles \n')
vel_prof_L = profile(:,4:6);
DCMnb_prof = profile(:,10:18);
delv_b = gendv(vel_prof_L,DCMnb_prof);
deltheta = gendthet(DCMnb_prof);

vel_L = profile(1,4:6);
pos(1,1:3) = profile(1,1:3);
dcmnb = [profile(1,10:12); profile(1,13:15); profile(1,16:18)];
dcmbn = dcmnb';

qua = dcm2qua(dcmbn);                % initialize quaternion

eroll(1) = roll(1); epit(1) = pitch(1); eyaw(1) = yaw(1);

fprintf(1,' Starting nav computations \n')
C = [0 1 0; 1 0 0; 0 0 -1];
for i = 2:npts,
   dtheta = deltheta(i-1,:);
   qua = quaupdat(qua,dtheta);        % perform quaternion updating
   DCMbn = qua2dcm(qua);              % convert quaternion to direction cosine matrix
     rpy = dcm2eulr(DCMbn);
     eroll(i) = rpy(1)*180/pi; 
     epit(i) = rpy(2)*180/pi;
     eyaw(i) = rpy(3)*180/pi;
   delv_L = C*DCMbn*delv_b(i-1,:)';
   vel_L(i,:) = vel_L(i-1,:) + delv_L';
   deltat = profile(i,19) - profile(i-1,19);
   pos(i,1:3)= pos(i-1,1:3) + deltat*mean([vel_L(i,1:3); vel_L(i-1,1:3)]);
end

plot3(pos(:,1),pos(:,2),pos(:,3))
axis equal
title('INSDEM04:  Quaternion Processing INS Computed Flight Path')
xlabel('east (meters)')
ylabel('north (meters)')
zlabel('up (meters)')
grid

figure
plot(time,pos(:,1)-profile(:,1),'o',...
   time,pos(:,2)-profile(:,2),'x',...
   time,pos(:,3)-profile(:,3),'+')
axis([0 80 -0.02 0.02])
title('Position Errors Due To Algorithm Limitations')
ylabel('error in meters')
xlabel('run time in seconds')
text(57,0.01,'east')
text(45,-0.01,'north')
text(62,-0.003,'up')
text(10,0.015,'Quaternion processing')
