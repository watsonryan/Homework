%  sndemo32.m      Compute Power Spectral Density of a C/A-code
%

clear all
close all

ca = prncode(10);        % Generate code for C/A PRN 10
ca2 = [];
for i = 1:max(size(ca)),
    ca2 = [ca2 ca(i) ca(i) ca(i) ca(i) ca(i)];    % Use a sample and hold to increase the
end                             % effective sampling frequency from 
                                % from 1.023 MHz to 5.115 MHz

[l,r] = gpscor(ca2,ca2,1022);
psdca = fft(r);

fs = 5*1.023e6;
num = max(size(r));
f = ((1:num)-1)*fs/num;
close
plot(f*0.001,20*log10(abs(psdca)+1))
axis([0 2.5*1.023e3 0 90])
title('C/A PRN 10 PSD')
ylabel('relative line component power in dB')
xlabel('frequency in kHz')

