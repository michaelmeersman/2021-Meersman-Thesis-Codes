clear; close all

width = 6;
height = 5;

% Pwelch parameters
fftsize = 2^15;
WINDOW_TYPE = hanning(fftsize);
NFFT = size(WINDOW_TYPE,1);
WINDOW_OVERLAP = fftsize/2;
SAMPLEFREQ = 50000;

% Apply band pass filter to signal
lowpass = 1e3;
highpass = 1e0;
d1 = designfilt('lowpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',lowpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
d2 = designfilt('highpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',highpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');

for throttle = [40 50 60 70 80 90]
    
    NI_data = readmatrix(['C:\Users\micha\Documents\' num2str(throttle) ' throttle.csv']);
    accel = NI_data(1:end,2);
    
    
    accel = filtfilt(d1,accel);
    accel = filtfilt(d2,accel);
    
    [power_mag, power_freq] = ...
        pwelch(accel, WINDOW_TYPE, WINDOW_OVERLAP, NFFT, SAMPLEFREQ);
    
    
    loglog(power_freq,power_mag,'DisplayName',[num2str(throttle) ' $\%$ Throttle'],'LineWidth',1.2)
    hold on
end

xlabel('Frequency, Hz','interpreter', 'latex')
ylabel('PSD','interpreter', 'latex')
set(gca, 'LineWidth',0.75 ,'FontSize',16)
grid on
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size
% xlim([highpass lowpass])
xlim([280 600])
% title('Wingtip Accelerometer Ground Tests','interpreter','latex')
legend('location','northeast','interpreter','latex','numcolumns',2)