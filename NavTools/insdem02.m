%
%   insdem02.m            Simplified 3-d navigation 
%
%     Vehicle starts from a stationary position, undergoes
%     constant acceleration for 10 seconds, then . . .
%
%     - Non-rotating, flat earth without gravity assumed.
%       Since the run time is extremely short (80 seconds),
%       the earth assumptions are valid.  'Without gravity'
%       is another way of saying that the gravity model is
%       perfect
%
%     - Locally-level, east-north-up coordinate system
%       (x: east, y: north, z: up)
%     - Heading (i.e., body attitude) is measured clockwise
%       from north
%
%     - Ideal, error-free measurements are being used in
%       order to focus on the INS algorithms

clear all
fprintf(1,' Creating flight profile \n')
insdem01
close all

npts = max(size(time));

fprintf(1,' Creating delta-v and delta-theta profiles \n')

vel_prof_L = profile(:,4:6);             % extracting velocity profile

DCMnb_prof = profile(:,10:18);           % extracting true nav-to-body
%                                        % direction cosine matrix

delv_b = gendv(vel_prof_L,DCMnb_prof);   % generating ideal, error-free
%                                        % delta-V measurements

deltheta = gendthet(DCMnb_prof);         % generating ideal, error-free
%                                        % delta-theta measurements
%                                        % (recall that the time span
%                                        % is very short so we can ignore
%                                        % the fact that we are actually
%                                        % moving over a rotating, 
%                                        % ellipsoidal earth

vel_L = profile(1,4:6);                % initializing velocity

pos(1,1:3) = profile(1,1:3);             % initializing position

%                                        % initializing vehicle attitude
DCMnb = [profile(1,10:12); profile(1,13:15); profile(1,16:18)];
DCMbn = DCMnb';

%                                        % setting first Euler angle
%                                        % measurements equal to true
eroll(1) = roll(1); epit(1) = pitch(1); eyaw(1) = yaw(1);
%estdcmbn(1,1:9) = [DCMbn(1,1:3) DCMbn(2,1:3) DCMbn(3,1:3)];
%trudcmbn(1,1:9) = [DCMbn(1,1:3) DCMbn(2,1:3) DCMbn(3,1:3)];

fprintf(1,' Starting nav computations \n')
C = [0 1 0; 1 0 0; 0 0 -1];       % Conversion between NED and ENU

for i = 2:npts,
   dtheta = deltheta(i-1,:);            % extract the current delta-theta
                                        % measurement from the profile
   
   DCMbn = bodupdat(DCMbn,dtheta);      % use the delta-theta to update
                                        % the body-to-nav direction
                                        % cosine matrix (DCM)
   
      rpy = dcm2eulr(DCMbn);            % Extract Euler angles from
      eroll(i) = rpy(1)*180/pi;         % the INS-derived DCM
      epit(i) = rpy(2)*180/pi;
      eyaw(i) = rpy(3)*180/pi;
      
   delv_L = C*DCMbn*delv_b(i-1,:)';     % Convert the delta-V measurement
                                        % from body coordinates to   
                                        % local-level coordinates
                                        
   vel_L(i,:) = vel_L(i-1,:) + delv_L'; % update local-level velocity
                                        
   deltat = profile(i,19) - profile(i-1,19);  % compute the time span of the
                                              % current measurement interval
                                              % (Note: simulation is multi-rate
                                              % so the intervals are not constant
   
                                % use traezoidal integration of velocity
                                % to update position
   pos(i,1:3)= pos(i-1,1:3) + deltat*mean([vel_L(i,1:3); vel_L(i-1,1:3)]);
end

plot3(pos(:,1),pos(:,2),pos(:,3))
axis equal
title('INSDEM02:  INS Computed Flight Path')
xlabel('east (meters)')
ylabel('north (meters)')
zlabel('up (meters)')
grid

figure
plot(time,pos(:,1)-profile(:,1),'o',...
   time,pos(:,2)-profile(:,2),'x',...
   time,pos(:,3)-profile(:,3),'+')
axis([0 80 -0.02 0.02])
title('INSDEM02: Position Errors Due To Algorithm Limitations')
ylabel('error in meters')
xlabel('run time in seconds')
text(57,0.01,'east')
text(45,-0.01,'north')
text(62,-0.003,'up')
