%
%   insdem21.m            Simplified 3-d navigation 
%
%     F-16 flight profile from INSDEM20 used here
%
%     - Non-rotating, flat earth without gravity assumed.
%       Since neither the run time nor the flight path are
%       particularly short (ten minutes and 70 kilometers),
%       the earth assumptions are not valid.  However, it is
%       easier to define a dynamic flight path first in the
%       local-level frame and then translate it to the earth
%       frame later.  See INSDEM27 for how this is done.
%
%       This example demonstrates the algorithmic errors since
%       there are no sensor errors being modeled.  
%
%       'Without gravity' is another way to say that the
%       gravity model is perfect
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
insdem20
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
eroll(1) = roll(1); epit(1) = pitch(1); eyaw(1) = yaw(1);
%estdcmbn(1,1:9) = [DCMbn(1,1:3) DCMbn(2,1:3) DCMbn(3,1:3)];
%trudcmbn(1,1:9) = [DCMbn(1,1:3) DCMbn(2,1:3) DCMbn(3,1:3)];

fprintf(1,' Starting nav computations \n')
C = [0 1 0; 1 0 0; 0 0 -1];       % Conversion between NED and ENU

h = waitbar(0,' Running through Time Loop');
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
   waitbar(i/npts,h)
end
close(h)

plot3(pos(:,1),pos(:,2),pos(:,3))
axis equal
title('INSDEM21:  INS Computed Flight Path')
xlabel('east (meters)')
ylabel('north (meters)')
zlabel('up (meters)')
grid

figure
plot(time,pos(:,1)-profile(:,1),'o',...
   time,pos(:,2)-profile(:,2),'x',...
   time,pos(:,3)-profile(:,3),'+')
title('INSDEM21: Position Errors Due To Algorithm Limitations')
ylabel('error in meters')
xlabel('run time in seconds')
legend('east error','north error','up error')
