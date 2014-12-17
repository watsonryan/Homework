% PRN cross-correlation example
%% cacode function from 
%% http://www.mathworks.com/matlabcentral/fileexchange/14670-gps-c-a-code-generator


% generate prn 1 C/A code
prn1=cacode(1);
% generate prn 2 C/A code
prn2=cacode(2);

% plot them just to provide a visual
subplot(311),stairs(prn1(1:100))
title('PRN Orthogonality Example')
legend('C/A Code PRN2')
ylim([-1 2])
subplot(312),stairs(prn2(1:100))
legend('C/A Code PRN2')
ylim([-1 2])

% compare the auto-correlation with PRN1 with the cross-correlation
% auto
Corr = (xcorr(prn1(1:end-9),prn1(10:end),'coeff'));
subplot(313),plot(xcorr(prn1(1:end-9),prn1(10:end),'coeff'),'r--','linewidth',2)
hold on
% cross
plot(xcorr(prn1(1:end-9),prn2(10:end),'coeff'),'linewidth',2)
ylim([0.4 1])
xlim([1023-100,1023+100])
legend('xcorr(PRN1,PRN1)','xcorr(PRN1,PRN2)')