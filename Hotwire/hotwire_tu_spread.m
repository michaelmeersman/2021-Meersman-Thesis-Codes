clear; close


for k = 1:5
    
    data(k) = load(['Flight_CVA_case_' num2str(k) '_adjusted.mat']);
    
end

%%

% figure
% subplot(5,1,1)
% plot(data(2).h_data(:,1))
% ylim([3.1 3.3])
% 
% subplot(5,1,2)
% plot(data(1).h_data(:,1))
% ylim([3.1 3.3])
% 
% subplot(5,1,3)
% plot(data(3).h_data(:,1))
% ylim([3.1 3.3])
% 
% subplot(5,1,4)
% plot(data(4).h_data(:,1))
% ylim([3.1 3.3])
% 
% subplot(5,1,5)
% plot(data(5).h_data(:,1))
% ylim([3.1 3.3])


% 
% for k = 1:5
%    
%     
%     velocity_mean(k) = mean(data(k).h_v);
%     
%     rms_velocity(k) = sqrt(1/length(data(k).h_v)*sum((data(k).h_v-velocity_mean(k)).^2));
%     
%     voltage_mean(k) = mean(data(k).h_data);
%     
%     rms_voltage(k) = sqrt(1/length(data(k).h_data)*sum((data(k).h_data-voltage_mean(k)).^2));
%     
%     Tu(k) = rms_velocity(k)/velocity_mean(k)*100;
% end
% 
% figure
% yyaxis left
% plot(rms_velocity)
% hold on
% yyaxis right
% plot(rms_voltage)
% 
% figure
% scatter(velocity_mean,Tu)

%% Split up data for each test sequence into 0.5s segments
close all

SAMPLEFREQ = 50000;
lowpass = 1e3;
highpass = 2e1;
d1 = designfilt('bandpassiir','FilterOrder',6, ...
    'HalfPowerFrequency1',highpass,'HalfPowerFrequency2',lowpass,...
    'SampleRate',SAMPLEFREQ,'DesignMethod','butter');

testcases = [1 2 3 4 5];
for k = testcases
    data(k).mov_av = movmean(data(k).h_v,SAMPLEFREQ);
    data(k).vel_bp = filtfilt(d1,data(k).h_v)+data(k).mov_av;
end


time_s = 0.2;
seglen = SAMPLEFREQ*time_s;

length_test = zeros(1,length(data));
nseg = zeros(1,length(data));

for k = testcases
    length_test(k) = floor(length(data(k).vel_bp)/seglen)*seglen;
    nseg(k) = length_test(k)/seglen;
end

for k = testcases
    data(k).mtx = zeros(nseg(k),seglen);
    data(k).seg_mean = zeros(nseg(k),1);
    
    for i = 1:nseg(k)
        data(k).mtx(i,:) = data(k).vel_bp((i-1)*seglen+1:i*seglen);
        data(k).seg_mean(i) = mean(data(k).mtx(i,:));
    end
end

for k = testcases
    data(k).seg_RMS = zeros(nseg(k),1);
    data(k).Tu = zeros(nseg(k),1);
    data(k).seg_prime = zeros(nseg(k),seglen);
    
    for i = 1:nseg(k)
        
        for q = 1:seglen
            data(k).seg_prime(i,q) = data(k).mtx(i,q) - data(k).seg_mean(i);
            data(k).prime_squared(i,q) = data(k).seg_prime(i,q)^2;
        end
        
        data(k).seg_RMS(i) = sqrt(1/seglen*sum(data(k).prime_squared(i,:)));
        
        data(k).rms_func(i) = rms(data(k).seg_prime(i,:));
        
        data(k).Tu(i) = data(k).rms_func(i)/data(k).seg_mean(i)*100; % Tu as percent for each speed
    end
    RMS_mean(k) = mean(data(k).seg_RMS(:));
    Tu_mean(k) = mean(data(k).Tu);
    Velocity_mean(k) = mean(data(k).seg_mean);
end

% Find the maximum and minimum Tu cases overall and generate energy spectra
max_case = 3;
min_case = 2;
[tu_max, max_loc] = max(data(max_case).Tu);
[tu_min, min_loc] = min(data(min_case).Tu);
RMS_overall = mean(RMS_mean(:));

fftsize = 2^13;
WINDOW_TYPE = hanning(fftsize);
WINDOW_OVERLAP = fftsize/2;

[power_mag, power_freq] = ...
        pwelch(data(max_case).mtx(max_loc,:), WINDOW_TYPE, WINDOW_OVERLAP, fftsize, SAMPLEFREQ);
max_power_mag = power_mag;
max_power_freq = power_freq;
[power_mag, power_freq] = ...
        pwelch(data(min_case).mtx(min_loc,:), WINDOW_TYPE, WINDOW_OVERLAP, fftsize, SAMPLEFREQ);
min_power_mag = power_mag;
min_power_freq = power_freq;

width = 5.25;
height = 4.5;

figure
loglog(max_power_freq,max_power_mag,'Linewidth',1.3)
grid on
hold on
loglog(min_power_freq,min_power_mag,'Linewidth',1.3)
x = linspace(1,10^4);
y = x.^(-5/3);
loglog(x,y,'k')
xlim([highpass-0.1*highpass lowpass])
xlab = xlabel('Frequency, Hz','interpreter', 'latex');
ylab = ylabel('PSD, $(m/s)^2/Hz$','interpreter', 'latex');
set(gca, 'LineWidth', 0.75 ,'FontSize',16)
% title(['CVA, 20 - 1000 Hz Band Pass, Tu = ' num2str(Tu,2) '\%'],'Interpreter', 'latex')
% title(['Tu = ' num2str(Tu,2) '\%'],'Interpreter', 'latex')
% title('CVA Hotwire','interpreter', 'latex')
% annotation('textbox', [0.72, 0.80, 0.1, 0.1], 'String', ['Tu = ' num2str(Tu,2) '$\%$'],'FontSize',20,'interpreter','latex','horizontalalignment','center','verticalalignment','top')
legend('MAX Tu','MIN Tu','$f^{-5/3}$ Line','Interpreter', 'latex')
xlab.FontSize = 20;
ylab.FontSize = 20;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size

figure
for k = testcases
    scatter(data(k).seg_mean,data(k).Tu,'Linewidth',1.2,'DisplayName',['Test Case ' num2str(testcases(k))])
    hold on
end
v_space = linspace(0,50);
scatter(Velocity_mean(testcases),Tu_mean(testcases),80,'k','Linewidth',1.7,'DisplayName','Averages')
avg_curve = plot(v_space,RMS_overall*100./v_space,'k','Linewidth',1.8);
avg_curve.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylim([0 0.9])
xlim([18 28.5])
ylab = ylabel('Tu, $\%$','Interpreter','latex');
xlab = xlabel('Velocity, m/s','Interpreter','latex');
legend('location','north','interpreter','latex','Fontsize',16,'numcolumns',3)
text(data(max_case).seg_mean(max_loc),data(max_case).Tu(max_loc)+0.0075,'\leftarrow MAX','Fontsize',20)
text(data(min_case).seg_mean(min_loc),data(min_case).Tu(min_loc)+0.0075,'\leftarrow MIN','Fontsize',20)
set(gca,'FontSize',16)
xlab.FontSize = 20;
ylab.FontSize = 20;
% title([num2str(time_s) ' second segments'],'Interpreter','latex','FontSize',20)

box on
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size


figure
choice = 2;
seconds = (1:length_test(choice))/50000;
plot(seconds,data(choice).h_v(1:length(seconds)))
hold on
plot(seconds,data(choice).vel_bp(1:length(seconds)),'Linewidth',1.2)
plot(seconds,data(choice).mov_av(1:length(seconds)),'Linewidth',1.5)
ylab = ylabel('Velocity, m/s','Interpreter','latex','FontSize',20);
xlab = xlabel('Time, s','Interpreter','latex','FontSize',20);
legend('Raw Hotwire Signal','Band-pass-filtered Signal',...
    'Moving Average of Signal','Interpreter','latex','FontSize',16)

title('$Re = 513k$','Interpreter','latex','FontSize',20)
% title([num2str(time_s) ' second segments'],'Interpreter','latex','FontSize',20)

% for i = 1:nseg(1)
%     xpos = seglen*i/50000;
%     xx = xline(xpos,'k','LineWidth',0.4);
%     xx.Annotation.LegendInformation.IconDisplayStyle = 'off';
% end
set(gca,'FontSize',16)
xlab.FontSize = 20;
ylab.FontSize = 20;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size







