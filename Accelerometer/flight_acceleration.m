clear; close

testcase = 5;
NI_data = readmatrix(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '1.1 Flight Test 9_18_21\CVA_Flight_3\CVA_' num2str(testcase) '.csv']);
accel = NI_data(100000:0.9*length(NI_data),2);


% Pwelch parameters
fftsize = 2^16;
WINDOW_TYPE = hanning(fftsize);
NFFT = size(WINDOW_TYPE,1);
WINDOW_OVERLAP = fftsize/2;
SAMPLEFREQ = 50000;


% Apply band pass filter to signal
lowpass = 1e3;
% highpass = 1e1;
d1 = designfilt('lowpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',lowpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
accel = filtfilt(d1,accel);


[power_mag, power_freq] = ...
        pwelch(accel, WINDOW_TYPE, WINDOW_OVERLAP, NFFT, SAMPLEFREQ);
    
loglog(power_freq,power_mag)
xlim([0 lowpass])
xlabel('Frequency, Hz','interpreter', 'latex')
ylabel('PSD','interpreter', 'latex')
title('Wingtip Accelerometer','interpreter', 'latex')
set(gca, 'FontName','Times','FontAngle','italic','FontWeight','bold')
set(gca, 'LineWidth',0.75 ,'FontSize',16)
grid on