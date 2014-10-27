%   Simplified 3-d navigation 

clear all

fprintf(1,' Creating flight profile \n')
%insdem32
load dem32_dat
close all

npts = max(size(time));

fprintf(1,' Creating delta-v and delta-theta profiles \n')

vel_prof_L = profile(:,4:6);             % extracting velocity profile

DCMnb_prof = profile(:,10:18);           % extracting true nav-to-body
%                                        % direction cosine matrix

fprintf(1,' Creating delta-v profile \n')
delv_b = gendv(vel_prof_L,DCMnb_prof);   % generating ideal, error-free
%                                        % delta-V measurements

fprintf(1,' Creating delta-theta profile \n')
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
eroll(1) = roll_deg(1); epit(1) = pitch_deg(1); eyaw(1) = yaw_deg(1);
%estdcmbn(1,1:9) = [DCMbn(1,1:3) DCMbn(2,1:3) DCMbn(3,1:3)];
%trudcmbn(1,1:9) = [DCMbn(1,1:3) DCMbn(2,1:3) DCMbn(3,1:3)];

fprintf(1,' Starting nav computations \n')
C = [0 1 0; 1 0 0; 0 0 -1];       % Conversion between NED and ENU

%
tru_vel = profile(:,4:6);
tru_vel_tim = profile(:,19);
%
h = waitbar(0,' Running through Time Loop');
for i = 2:npts,
    %
   dtheta = deltheta(i-1,:);            % extract the current delta-theta
                                        % measurement from the profile
   
   DCMbn = bodupdat(DCMbn,dtheta);      % use the delta-theta to update
                                        % the body-to-nav direction
                                        % cosine matrix (DCM)
   
      rpy = dcm2eulr(DCMbn);            % Extract Euler angles from
      eroll(i) = rpy(1)*180/pi;         % the INS-derived DCM
      epit(i) = rpy(2)*180/pi;
      eyaw(i) = rpy(3)*180/pi;
      
      if eyaw(i) > 180
          eyaw(i)=eyaw(i)-180;
      end
      
   delv_L = C*DCMbn*delv_b(i-1,:)';     % Convert the delta-V measurement
                                        % from body coordinates to   
                                        % local-level coordinates
                                        
   vel_L(i,:) = vel_L(i-1,:) + delv_L'; % update local-level velocity
                                        
   deltat = profile(i,19) - profile(i-1,19);  % compute the time span of the
                                              % current measurement interval
                                              % (Note: simulation is multi-rate
                                              % so the intervals are not constant

%                                 % use trapezoidal integration of velocity
%                                 % to update position
    pos(i,1:3)= pos(i-1,1:3) + deltat*mean([vel_L(i,1:3); vel_L(i-1,1:3)]);
    
    if yaw_deg(i) > 180
        yaw_deg(i)=yaw_deg(i)-360;
    end
%
   waitbar(i/npts,h)
end
close(h)


% figure
% plot3(profile(:,1),profile(:,2),profile(:,3),'r')
% axis equal
% title('Flight Path')
% xlabel('east (meters)')
% ylabel('north (meters)')
% zlabel('up (meters)')
% grid

figure
subplot(311)
plot(time,profile(:,1))
title('INSDEM32: Position Components')
ylabel('east in meters')
subplot(312)
plot(time,profile(:,2))
ylabel('north in meters')
subplot(313)
plot(time,profile(:,3))
ylabel('vertical in meters')
xlabel('run time in seconds')

figure
subplot(311)
plot(time,profile(:,4))
title('INSDEM32: Velocity Components')
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
plot(time_eulr,roll_deg)
title('INSDEM32: Euler Angles')
ylabel('roll angle in deg')
subplot(312)
plot(time_eulr,pitch_deg)
ylabel('pitch angle in deg')
subplot(313)
plot(time_eulr,yaw_deg)
ylabel('yaw angle in deg')
xlabel('run time in seconds')





