clear; close all

width = 5;
height = 5;
testcase = 5;
% Select which case to evaluate
flow = 10;
h_vel = importdata(['CVA_' num2str(testcase) '_adjusted.csv']);

% Import Pressure data for comparisson to pitot
flight_no = 3;
daq = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '1.1 Flight Test 9_18_21\RPI_Flight_' num2str(flight_no) '\pressure_data_case_' num2str(testcase) '_f.csv']);
pitot = daq(:,2);

% Parameters for calcTurbValues function
fftsize = 2^14;
WINDOW_TYPE = hanning(fftsize);
WINDOW_OVERLAP = fftsize/2;
SAMPLEFREQ = 50000;
NFFT = size(WINDOW_TYPE,1);


h_vel = mean(pitot)/mean(h_vel).*h_vel;
h_trace = h_vel;
h_mean = mean(h_vel);

% Apply band pass filter to signal
lowpass = 1e3;
highpass = 2e1;
d1 = designfilt('lowpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',lowpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
d2 = designfilt('highpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',highpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
h_vel = filtfilt(d1,h_vel-mean(h_vel)) + h_mean;
h_vel = filtfilt(d2,h_vel-mean(h_vel)) + h_mean;

% May need to cut off the first and last chunk of data to avoid distortion
% from "0-padding" due to ends of data set. Higher filter order, larger



h_prime = zeros(length(h_vel),1);
for j = 1:length(h_vel)
    h_prime(j) = (h_vel(j) - h_mean);
end
h_RMS = sqrt(1/length(h_vel)*sum(h_prime.^2));
Tu = h_RMS/h_mean*100; % Tu as percent for each speed


[wt,f] = cwt(h_prime,SAMPLEFREQ,'voicesperocatave',20);
wt = abs(wt);
t = (0:length(h_prime)-1)/SAMPLEFREQ;

[Tq,Fq] = meshgrid(0:14/1000:14,10:2000/length(f):2000);
WTq = interp2(t,f,wt,Tq,Fq);

pcolor(Tq,Fq,WTq)
ylim([20 1000])
xlabel('Time, s','interpreter','latex')
ylabel('Frequency, Hz','interpreter','latex')
c = colorbar('eastoutside');
c.Label.String = 'Magnitude';
c.Label.Interpreter = 'latex';
set(gca,'Fontsize',16,'clim',[0.01 0.04])
shading flat
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*130, height*120]); %<- Set size




