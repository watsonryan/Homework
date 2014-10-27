%
%  plot_kf_fb.m
%    plot kf-with-feedback example
%
clear all
close all
load kf_fb_out

figure
subplot(311)
plot(time/60,east_pos_err_KF,'b-',...
    time/60,x_rms_KF,'r--',time/60,-x_rms_KF,'r--')
axis([0 time_end -2 2])
title('East/North/Up Position Error in Meters')
ylabel('east')

subplot(312)
plot(time/60,north_pos_err_KF,'b-',...
    time/60,y_rms_KF,'r--',time/60,-y_rms_KF,'r--')
axis([0 time_end -2 2])
ylabel('north')

subplot(313)
plot(time/60,up_pos_err_KF,'b-',...
    time/60,z_rms_KF,'r--',time/60,-z_rms_KF,'r--')
axis([0 time_end -2 2])
ylabel('up')
xlabel('time in minutes')

figure
subplot(311)
plot(time/60,vel_l_err_KF(:,1),'b-',...
    time/60,x_v_rms_KF,'r--',time/60,-x_v_rms_KF,'r--')
axis([0 time_end -0.1 0.1])
title('Velocity Error in Meters/Second')
ylabel('east')

subplot(312)
plot(time/60,vel_l_err_KF(:,2),'b-',...
    time/60,y_v_rms_KF,'r--',time/60,-y_v_rms_KF,'r--')
axis([0 time_end -0.1 0.1])
ylabel('north')

subplot(313)
plot(time/60,vel_l_err_KF(:,3),'b-',...
    time/60,z_v_rms_KF,'r--',time/60,-z_v_rms_KF,'r--')
axis([0 time_end -0.1 0.1])
ylabel('up')
xlabel('time in minutes')

figure
subplot(311)
plot(time/60,roll_err_KF*(180/pi),'b-',...
    time/60,x_psi_rms_KF*(180/pi),'r--',time/60,-x_psi_rms_KF*(180/pi),'r--')
axis([0 time_end -.001 .001])
title('Attitude/PSI Estimation Error in Degrees')
ylabel('x-comp / roll')

subplot(312)
plot(time/60,pitch_err_KF*(180/pi),'b-',...
    time/60,y_psi_rms_KF*(180/pi),'r--',time/60,-y_psi_rms_KF*(180/pi),'r--')
ylabel('y-comp / pitch')
axis([0 time_end -.001 .001])

subplot(313)
plot(time/60,yaw_err_KF*(180/pi),'b-',...
    time/60,z_psi_rms_KF*(180/pi),'r--',time/60,-z_psi_rms_KF*(180/pi),'r--')
ylabel('z-comp / yaw')
xlabel('time in minutes')
axis([0 time_end -.01 .01])

g = gravity(0,0);
mpss2microg = 1/(g*1e-6);  % meters-per-seconds-squared to micro-g
figure
subplot(311)   % Note: truth is given in vxbias in m/s^2
plot(time/60,(accum_accel_x+state_est(10,:)'-vxbias)*mpss2microg,'b-',...
    time/60,x_accel_rms_KF*mpss2microg,'r--',...
    time/60,-x_accel_rms_KF*mpss2microg,'r--')
axis([0 time_end -15 15])
title('Accel Bias Estimation Error in Micro-g')
ylabel('x accel')

subplot(312)
plot(time/60,(accum_accel_y+state_est(11,:)'-vybias)*mpss2microg,'b-',...
    time/60,y_accel_rms_KF*mpss2microg,'r--',...
    time/60,-y_accel_rms_KF*mpss2microg,'r--')
axis([0 time_end -15 15])
ylabel('y accel')

subplot(313)
plot(time/60,(accum_accel_z+state_est(12,:)'-vzbias)*mpss2microg,'b-',...
    time/60,z_accel_rms_KF*mpss2microg,'r--',...
    time/60,-z_accel_rms_KF*mpss2microg,'r--')
axis([0 time_end -15 15])
ylabel('z accel')
xlabel('time in minutes')

rps2dph = (180/pi)*3600;  % conversion constant from rad/sec to deg/hr
figure
subplot(311)   % NOTE: Truth is given in thxbias in deg/hr
plot(time/60,(accum_gyro_x+state_est(13,:)')*rps2dph-thxbias,'b-',...
    time/60,x_gyro_rms_KF*rps2dph,'r--',...
    time/60,-x_gyro_rms_KF*rps2dph,'r--')
title('Gyro Bias Estimation Error in Deg/Hr')
ylabel('x gyro')

subplot(312)
plot(time/60,(accum_gyro_y+state_est(14,:)')*rps2dph-thybias,'b-',...
    time/60,y_gyro_rms_KF*rps2dph,'r--',...
    time/60,-y_gyro_rms_KF*rps2dph,'r--')
ylabel('y gyro')

subplot(313)
plot(time/60,(accum_gyro_z+state_est(15,:)')*rps2dph-thzbias,'b-',...
    time/60,z_gyro_rms_KF*rps2dph,'r--',...
    time/60,-z_gyro_rms_KF*rps2dph,'r--')
ylabel('z gyro')
xlabel('time in minutes')

figure
subplot(311)   % NOTE: Truth is given in thxbias in deg/hr
plot(time/60,(accum_gyro_x+state_est(13,:)')*rps2dph-thxbias,'b-',...
    time/60,x_gyro_rms_KF*rps2dph,'r--',...
    time/60,-x_gyro_rms_KF*rps2dph,'r--')
axis([0 time_end -0.02 0.02])
title('Gyro Bias Estimation Error in Deg/Hr')
ylabel('x gyro')

subplot(312)
plot(time/60,(accum_gyro_y+state_est(14,:)')*rps2dph-thybias,'b-',...
    time/60,y_gyro_rms_KF*rps2dph,'r--',...
    time/60,-y_gyro_rms_KF*rps2dph,'r--')
axis([0 time_end -0.02 0.02])
ylabel('y gyro')

subplot(313)
plot(time/60,(accum_gyro_z+state_est(15,:)')*rps2dph-thzbias,'b-',...
    time/60,z_gyro_rms_KF*rps2dph,'r--',...
    time/60,-z_gyro_rms_KF*rps2dph,'r--')
axis([0 time_end -0.02 0.02])
ylabel('z gyro')
xlabel('time in minutes')

