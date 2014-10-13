%% INS Simulation

run('FlightPath')

g=9.81;
C_update = [6.12323399573677e-17,-1,0;1,6.12323399573677e-17,0;0,0,1];
V_update = [0.0375000000000000,2.29621274840129e-18,0];
P_update = [100.000187500000,200,50];
x(:,1) = 0; y(:,1) = 0; z(:,1) = 1.5708;

ENU=[0 1 0; 1 0 0; 0 0 -1]; 

h = waitbar(0,' Running through Time Loop');
for i = 2:npts,
    C_old = C_update;
    V_old = V_update(i-1,:);
    P_old = P_update(i-1,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%  Computes Skew-symmetric Matrix %%%%%%%%%%%%%%%%%%
    
    dtheta = deltheta(i-1,:);
    sigma_x = dtheta(1);
    sigma_y = dtheta(2);
    sigma_z = dtheta(3);
    
    S(:,:,i) = [  0      -sigma_z  sigma_y;
                sigma_z    0       -sigma_x;
                -sigma_y  sigma_x    0     ];
            
    
    magn = norm(dtheta);
    
    if magn == 0,
     Sbb = eye(3);
    else
     Sbb = eye(3) + (sin(magn)/magn)*S(:,:,i) + ( (1-cos(magn))/magn^2 )*S(:,:,i)*S(:,:,i);
    end
           
            
    C_update = C_old*Sbb;
    C_avg = .5*(C_update + C_old);
    delv=ENU*C_avg*delv_b(i-1,1:3)';
    V_update(i,1:3) = V_old + delv';
    P_update(i,1:3) = P_old + (((V_old +V_update(i,1:3)))*(deltat/2));
    
    
     x(:,i) = (atan2(C_update(3,2),C_update(3,3)));  
     y(:,i) = (asin(-C_update(3,1)));
     z(:,i) = (atan2(C_update(2,1),C_update(1,1)));
    
  

    waitbar(i/npts,h)
end
    close(h)
    
    
    x=x*(180/pi);
    y=y*(180/pi);
    z=z*(180/pi);
    
    
% figure
% plot3(P_update(:,1),P_update(:,2),P_update(:,3),'r')
% axis equal
% title('Position Solution')
% xlabel('x (meters)')
% ylabel('y (meters)')
% zlabel('z (meters)')
% grid
    
    
figure
subplot(311)
plot(time_eulr,x)
title('Euler Angles Inertial')
ylabel('roll angle in deg')
subplot(312)
plot(time_eulr,y)
ylabel('pitch angle in deg')
subplot(313)
plot(time_eulr,z)
ylabel('yaw angle in deg')
xlabel('run time in seconds')
