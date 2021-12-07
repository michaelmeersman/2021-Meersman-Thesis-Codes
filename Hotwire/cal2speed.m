%% Load all of the calibration data
clear; close

spd = 11:30;
avg = zeros(length(spd),1);

q = 1;
for i = spd
    data = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\' ...
        '2. DAQ\4. Tao Anemometer\Hotwire Calibration 07_07_21\07_07_21_P01\' ...
        '07_07_21 Calibration Data\CVA_cal_' num2str(i) '.csv'],'\t');
    cal_data = data.data;
    avg(q) = mean(cal_data(:,1));
    q = q + 1;
end
cases = readmatrix(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\' ...
    '2. DAQ\4. Tao Anemometer\Hotwire Calibration 07_07_21\' ...
    '07_07_21_P01\07_07_21_cases.xlsx'],'Range','C2:C21');
    
%% Create the polynomial fit and check for accuracy
p = polyfit(avg,cases,4);

scatter(avg,cases,'k')
hold on
range = 3:0.01:3.3;
fit = polyval(p,range);
plot(range,fit,'r')
xlim([3 3.3])
ylim([10 30])
set(gca, 'FontName','Times New Roman','FontSize',14)
ylabel('Calibrator speed, m/s','interpreter','latex','fontsize',20)
xlabel('CVA Output Voltage','interpreter','latex','fontsize',20)
legend('Hotwire Average','4th-Order Fit','Location','northwest','interpreter','latex','fontsize',16)
grid on
hold off

%% Apply polynomial fit to a dataset

testcase = 5;
% Create a new .csv file to store the calibrated data
f_id = fopen(['Flight_CVA_case_' num2str(testcase) '_adjusted.csv'],'w+');
% Open the uncalibrated data set for the selected flow case
h_data = readmatrix(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '1.1 Flight Test 9_18_21\CVA_Flight_3\CVA_' num2str(testcase) '.csv']);
% h_data = h_data.data;

% Locations in data where trigger signal vanishes - 0.5s of data
region = [553522;
          651420;
          925131;
          446285;
          913579]-25000;
      
h_data = h_data(25000:region(testcase),1);% Start 0.5s after trigger signal appears

% Apply the voltage to flow speed conversion using the polynomial values
% obtained from calibration
h_v = zeros(length(h_data),1);
for i = 1:length(h_data) 
   h_v(i) =  p(1)*h_data(i)^4 + p(2)*h_data(i)^3 + p(3)*h_data(i)^2 ...
       + p(4)*h_data(i) + p(5);
end
fprintf(f_id,'%2.8f\r\n',h_v);
fclose(f_id);

save(['Flight_CVA_case_' num2str(testcase) '_adjusted.mat'])

%% Compute Turbulent Values

% Select which case to evaluate
flow = 10;
h_vel = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\' ...
    '2. DAQ\4. Tao Anemometer\CVA_' num2str(flow) ...
    '_m_s_adjusted.csv'],'\t');

% Parameters for calcTurbValues function
fftsize = 2^18;
% Higher resolution for larger fftsize. Increasing resolution here leads to
% more noise as it reduces the number of windows you can obtain. Number of
% windows is equal to length of sample size divided by fftsize times 2
% rounded to next complete window size. (Can only do whole windows).
% "Accumulates" nearby frequency content into the same peak -- disadvantage
% is that you may see a peak that is not accurate. 
WINDOW_TYPE = hanning(fftsize);
WINDOW_OVERLAP = fftsize/2;
% Reducing/increasing number of windows to consider with pwelch -- larger
% number of windows decreases noise (above 10 is good rule of thumb)
% Having an overlap of 50% yields the maximum number of windows I can get
% without duplicating information -- not actually generating real data
SAMPLEFREQ = 100000;

% Subtract mean flow from velocity readings to obtain fluctuations
h_fluct = h_vel;

% Apply band pass filter to signal
lowpass = 1e3;
highpass = 2e1;
d1 = designfilt('lowpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',lowpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
d2 = designfilt('highpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',highpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
h_vel = filtfilt(d1,h_vel-mean(h_vel))+ mean(h_vel);
h_vel = filtfilt(d2,h_vel-mean(h_vel))+ mean(h_vel);

% May need to cut off the first and last chunk of data to avoid distortion
% from "0-padding" due to ends of data set. Higher filter order, larger
% error may occur on beginning and end windows
h_vel = h_vel(100:length(h_vel)-100);

[Tu, Vel_RMS, power_mag, power_freq] = calcTurbValues(1.65E-5,h_vel,SAMPLEFREQ,WINDOW_TYPE,WINDOW_OVERLAP,fftsize);

loglog(power_freq,power_mag)
grid on
xlim([highpass-0.1*highpass lowpass])

