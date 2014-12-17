%% INS Simulation
%%Authors: Ryan Watson

%%Last Edited:   10/30/2014

%run('FlightPath')

C_update = [6.12323399573677e-17,-1,0;1,6.12323399573677e-17,0;0,0,1];
V_update = [0.0212132034355964,0.0212132034355964,0];
P_update = [0.000106066017177982,0.000106066017177982,0];
x(:,1) = 0; y(:,1) = 0; z(:,1) = 1.5708;

ENU=[0 1 0; 1 0 0; 0 0 -1]; 

h = waitbar(0,' Running through Time Loop');
for i = 2:length(time)
    C_old = C_update;
    V_old = V_update(i-1,:);
    P_old = P_update(i-1,:);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%  Computes Skew-symmetric Matrix %%%%%%%%%%%%%%%%%%
    
    dtheta = est_dtheta(i-1,:);
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
    delv=ENU*C_update*est_dv(i-1,1:3)';
    V_update(i,1:3) = V_old + (deltav_b(i,:));
    P_update(i,1:3) = P_old + (((V_old + V_update(i,1:3)))*(.005));
   
    
    
     x(:,i) = (atan2(C_update(3,2),C_update(3,3)));  
     y(:,i) = (asin(-C_update(3,1)));
     z(:,i) = (atan2(C_update(2,1),C_update(1,1)));
    
  

    waitbar(i/npts,h)
end
    close(h)
    
    
    x=x*(180/pi);
    y=y*(180/pi);
    z=z*(180/pi);
    
    
    for i=1:length(yaw_deg);
    if yaw_deg(i)>180
        yaw_deg(i)=yaw_deg(i)-360;
    else
        yaw_deg(i)=yaw_deg(i);
    end
    end
    

figure
subplot(311)
plot(time,P_update(:,1)-profile(:,1))
title('Position Components Inertial')
ylabel('east in meters')
subplot(312)
plot(time,P_update(:,2)-profile(:,2))
ylabel('north in meters')
subplot(313)
plot(time,P_update(:,3)-profile(:,3))
ylabel('vertical in meters')
xlabel('run time in seconds')

figure
subplot(311)
plot(time,V_update(:,1)-profile(:,4))
title(' Velocity Components Inertial')
ylabel('east vel in m/s')
subplot(312)
plot(time,V_update(:,2)-profile(:,5))
ylabel('north vel in m/s')
subplot(313)
plot(time,V_update(:,3)-profile(:,6))
ylabel('vertical vel in m/s')
xlabel('run time in seconds')   

figure
subplot(311)
plot(time_eulr,x-roll_deg)
title('Euler Angles Inertial')
ylabel('roll angle in deg')
subplot(312)
plot(time_eulr,y-pitch_deg)
ylabel('pitch angle in deg')
subplot(313)
plot(time_eulr,z-yaw_deg)
ylabel('yaw angle in deg')
xlabel('run time in seconds')
