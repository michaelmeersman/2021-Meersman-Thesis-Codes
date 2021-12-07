clear; close

width = 6;
height = 5;

NI_data = readmatrix(['C:\Users\micha\Documents\Wing Acceleration Ground Tests\bending fan off.csv']);
accel = NI_data(1:end,2);

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


accel = filtfilt(d1,accel);
accel = filtfilt(d2,accel);

[power_mag, power_freq] = ...
    pwelch(accel, WINDOW_TYPE, WINDOW_OVERLAP, NFFT, SAMPLEFREQ);

[maxx, loc] = max(power_mag);


loglog(power_freq,power_mag,'DisplayName','Wingtip Accelerometer','LineWidth',1.2)
hold on

xlabel('Frequency, Hz','interpreter', 'latex')
ylabel('PSD','interpreter', 'latex')
set(gca, 'LineWidth',0.75 ,'FontSize',16)
grid on
xxxx = xline(power_freq(loc),'Label',['Wing Eigenfrequency ' num2str(power_freq(loc),3) ' Hz'],'LabelHorizontalAlignment','center','LabelVerticalAlignment','middle','LineWidth',1.5,'FontSize',14);
xxxx.Annotation.LegendInformation.IconDisplayStyle = 'off';
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size
xlim([highpass lowpass])
% title('Bending Test Frequency Spectrum','interpreter','latex')
legend('location','northeast','interpreter','latex')