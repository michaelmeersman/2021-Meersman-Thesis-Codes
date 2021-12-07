clear; close all

testcase = 2;
% Select which case to evaluate
h_vel = importdata(['CVA_' num2str(testcase) '_adjusted.csv']);

% Parameters for calcTurbValues function
fftsize = 2^16;
WINDOW_TYPE = hanning(fftsize);
WINDOW_OVERLAP = fftsize/2;
SAMPLEFREQ = 50000;

Re = [0 513 395 0 378];

h_mean = mean(h_vel);

% Apply band pass filter to signal
lowpass = 1e3;
highpass = 2e1;

d1 = designfilt('lowpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',lowpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
d2 = designfilt('highpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',highpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
h_vel = filtfilt(d1,h_vel-mean(h_vel)) + mean(h_vel);
h_vel = filtfilt(d2,h_vel-mean(h_vel)) + mean(h_vel);

[Tu, Vel_RMS, power_mag, power_freq] = calcTurbValues(1.65E-5,h_vel,SAMPLEFREQ,WINDOW_TYPE,WINDOW_OVERLAP,fftsize);

h_prime = zeros(length(h_vel),1);
for j = 1:length(h_vel)
    h_prime(j) = (h_vel(j) - h_mean);
end
h_RMS = sqrt(1/length(h_vel)*sum(h_prime.^2));
Tu = h_RMS/h_mean*100; % Tu as percent for each speed

width = 6;
height = 5;

figure
loglog(power_freq,power_mag)
grid on
hold on
% xln = xline(436,'-.r','Label','436 Hz','interpreter', 'latex','FontSize',20);
% xln.Annotation.LegendInformation.IconDisplayStyle = 'off';
x = linspace(1,10^4);
y = x.^(-5/3);
loglog(x,y,'k')
xlim([highpass-0.1*highpass lowpass])
xlabel('Frequency, Hz','interpreter', 'latex')
ylabel('PSD, $(m/s)^2/Hz$','interpreter', 'latex')
set(gca, 'LineWidth', 0.75 ,'FontSize',20)
% title(['CVA, 20 - 1000 Hz Band Pass, Tu = ' num2str(Tu,2) '\%'],'Interpreter', 'latex')
% title(['Tu = ' num2str(Tu,2) '\%'],'Interpreter', 'latex')
title(['$Re = $' num2str(Re(testcase)) 'k'],'interpreter', 'latex')
% annotation('textbox', [0.72, 0.80, 0.1, 0.1], 'String', ['Tu = ' num2str(Tu,2) '$\%$'],'FontSize',20,'interpreter','latex','horizontalalignment','center','verticalalignment','top')
legend('Hotwire Signal','$f^{-5/3}$ Line','location','northwest','Interpreter', 'latex')
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size




