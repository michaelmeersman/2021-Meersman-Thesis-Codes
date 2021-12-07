clear; close all

testcase = 2;
% Select which case to evaluate
flow = 10;
h_vel = importdata(['CVA_' num2str(testcase) '_adjusted.csv']);
% h_vel = importdata(['Flight_CVA_case_' num2str(testcase) '_adjusted.csv']);

% Import Pressure data for comparisson to pitot
flight_no = 3;
daq = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '1.1 Flight Test 9_18_21\RPI_Flight_' num2str(flight_no) '\pressure_data_case_' num2str(testcase) '_f.csv']);
pitot = daq(:,2);
alpha = round(mean(daq(:,1))-3,1);
v_ind = mean(pitot);            % Average indicated airspeed, m/s
mu = 1.894e-5;                  % Dynamic viscosity, m/s^2
rho = 1.13;                     % Air density, kg/m^3
chord = 0.305;                      % Chord length, m
RE = round(rho*v_ind*cosd(35)*chord/mu,-3);   % Reynolds number based on indicated 
                                % airspeed in chordwise direction

% Parameters for calcTurbValues function
fftsize = 2^16;
WINDOW_TYPE = hanning(fftsize);
WINDOW_OVERLAP = fftsize/2;
SAMPLEFREQ = 50000;
numWINDOWS = length(h_vel)/fftsize;


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


[Tu, Vel_RMS, power_mag, power_freq] = calcTurbValues(1.65E-5,h_vel,SAMPLEFREQ,WINDOW_TYPE,WINDOW_OVERLAP,fftsize);

h_prime = zeros(length(h_vel),1);
for j = 1:length(h_vel)
    h_prime(j) = (h_vel(j) - h_mean);
end
h_RMS = sqrt(1/length(h_vel)*sum(h_prime.^2));
Tu = h_RMS/h_mean*100; % Tu as percent for each speed

width = 12;
height = 3.8;

figure
subplot(1,3,1)
loglog(power_freq,power_mag)
grid on
hold on
x = linspace(1,10^4);
y = x.^(-5/3);
loglog(x,y,'k')
xlim([highpass-0.1*highpass lowpass])
xlabel('Frequency, Hz','interpreter', 'latex')
ylabel('PSD','interpreter', 'latex')
xln = xline(435,'-.r','Linewidth',0.5);
xln.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(gca, 'LineWidth', 0.75 ,'FontSize',16)
% title(['CVA, 20 - 1000 Hz Band Pass, Tu = ' num2str(Tu,2) '\%'],'Interpreter', 'latex')
% title(['Tu = ' num2str(Tu,2) '\%'],'Interpreter', 'latex')
title('Hotwire','interpreter', 'latex')
% annotation('textbox', [0.72, 0.80, 0.1, 0.1], 'String', ['Tu = ' num2str(Tu,2) '$\%$'],'FontSize',20,'interpreter','latex','horizontalalignment','center','verticalalignment','top')
legend('Hotwire Signal','$f^{-5/3}$ Line','Interpreter', 'latex')
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*130, height*120]); %<- Set size

save(['CVA_case_' num2str(testcase) '.mat'])

NI_data = readmatrix(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '1.1 Flight Test 9_18_21\CVA_Flight_3\CVA_' num2str(testcase) '.csv']);
accel = NI_data(length(NI_data)*0.1:length(NI_data)*0.9,2);


% Pwelch parameters
fftsize = 2^16;
WINDOW_TYPE = hanning(fftsize);
NFFT = size(WINDOW_TYPE,1);
WINDOW_OVERLAP = fftsize/2;
SAMPLEFREQ = 50000;


% Apply band pass filter to signal
accel = filtfilt(d1,accel);
accel = filtfilt(d2,accel);

[power_mag, power_freq] = ...
        pwelch(accel, WINDOW_TYPE, WINDOW_OVERLAP, NFFT, SAMPLEFREQ);

subplot(1,3,2)
loglog(power_freq,power_mag)
xlim([0 lowpass])
xln = xline(434,'-.r','Linewidth',0.5);
xln.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlabel('Frequency, Hz','interpreter', 'latex')
ylabel('PSD','interpreter', 'latex')
title('Wingtip Accelerometer','Interpreter', 'latex')
set(gca, 'LineWidth',0.75 ,'FontSize',16)
grid on
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*130, height*120]); %<- Set size
xlim([highpass lowpass])



% Apply crossPSD
fftsize_C = 2^15;
WINDOW_TYPE_C = hanning(fftsize_C);
NFFT_C = size(WINDOW_TYPE_C,1);
WINDOW_OVERLAP_C = fftsize_C/2;

[cxy,cfr] = cpsd(accel,h_vel,WINDOW_TYPE_C, WINDOW_OVERLAP_C, NFFT_C, SAMPLEFREQ);

subplot(1,3,3)
loglog(cfr,abs(cxy))
xlim([highpass lowpass])
ylim([1e-10 1e-5])
xln = xline(434,'-.r','Linewidth',0.5);
xln.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlabel('Frequency, Hz','interpreter', 'latex')
ylabel('Cross-PSD','interpreter', 'latex')
title('Cross-PSD','Interpreter', 'latex')
set(gca, 'LineWidth',0.75 ,'FontSize',16)
grid on

% sgtitle(['Re = ' num2str(RE) ', $\alpha = ' num2str(alpha) '^{\circ}$'],'interpreter','latex','FontSize',20) 

%% Split hotwire signal into 20 segments and plot Tu vs vel

close

nseg = 20;
seglen = round(length(h_vel)/nseg)-1;

mtx = zeros(nseg,seglen);
seg_mean = zeros(nseg,1);
seg_RMS = zeros(nseg,1);
Tu = zeros(nseg,1);
seg_prime = zeros(nseg,seglen);

for i = 1:nseg
    
    mtx(i,:) = h_vel((i-1)*seglen+1:i*seglen);
    seg_mean(i) = mean(mtx(i,:));
    for k = 1:seglen
        seg_prime(i,k) = mtx(i,k) - seg_mean(i);
    end
    
    seg_RMS(i) = sqrt(1/length(seg_prime(i))*sum(seg_prime(i).^2));
    Tu(i) = seg_RMS(i)/seg_mean(i)*100; % Tu as percent for each speed
end

seg_time = round(seglen/50000,2);

figure
scatter(seg_mean,Tu)
xlabel('Speed, m/s','interpreter', 'latex')
ylabel('Tu, \%','interpreter', 'latex')
xlim([min(seg_mean) max(seg_mean)])
title(['Tu vs Speed for ' num2str(seg_time) ' s data segments'],'interpreter', 'latex')
set(gca, 'LineWidth',0.75 ,'FontSize',16)

%% Compare velocity traces of pitot and hotwire
close
sample_time = floor(length(h_vel)/SAMPLEFREQ);
pitot_clipped = pitot(length(pitot)*0.1:length(pitot)*0.9);
x_pitot = (0:1/(length(pitot_clipped)-1):1)*sample_time;
h_trace_clipped = h_trace(1:length(h_trace));
x_hotwire = 0:1/length(h_trace_clipped):1;
x_hotwire = (x_hotwire(1:end-1))*sample_time;
accel_clipped = accel(1:length(h_trace));

figure
subplot(2,1,1)
plot(x_hotwire,h_trace_clipped)
hold on
plot(x_pitot,pitot_clipped)
legend('CVA','Pitot','interpreter','latex')
% xlabel('Time, $s$','interpreter','latex')
ylabel('Velocity, $m/s$','interpreter','latex')
set(gca, 'FontSize', 16)
grid on
subplot(2,1,2)
plot(x_hotwire,accel_clipped)
pos = get(gcf, 'Position');
legend('Wingtip Accelerometer','interpreter','latex')
xlabel('Time, $s$','interpreter','latex')
ylabel('Acceleration, $m/s^2$','interpreter','latex')
set(gca, 'FontSize', 16)
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*80, height*140]); %<- Set size


sgtitle(['$Re_{mean}$ = ' num2str(RE) ', $\alpha_{mean} = ' num2str(alpha) '^{\circ}$'],'interpreter','latex','FontSize',20) 


%% Create spectrogram for hotwire signal

type = accel;
[wt,f] = cwt(type,SAMPLEFREQ,'voicesperocatave',12);
wt = abs(wt);
t = (0:length(type)-1)/SAMPLEFREQ;

[Tq,Fq] = meshgrid(0:floor(length(type)/50000)/1000:floor(length(type)/50000),10:2000/length(f):2000);
WTq = interp2(t,f,wt,Tq,Fq);

figure
pcolor(Tq,Fq,WTq)
ylim([20 1000])
xlabel('Time, s','interpreter','latex')
ylabel('Frequency, Hz','interpreter','latex')
c = colorbar('eastoutside');
c.Label.String = 'Magnitude';
c.Label.Interpreter = 'latex';
set(gca,'Fontsize',16)
shading flat
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*80, height*120]); %<- Set size

sgtitle(['Re = ' num2str(RE) ', $\alpha = ' num2str(alpha) '^{\circ}$'],'interpreter','latex','FontSize',20) 



