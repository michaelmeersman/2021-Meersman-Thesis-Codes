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



for k = 1:5
   
    
    velocity_mean(k) = mean(data(k).h_v);
    
    rms_velocity(k) = sqrt(1/length(data(k).h_v)*sum((data(k).h_v-velocity_mean(k)).^2));
    
    voltage_mean(k) = mean(data(k).h_data);
    
    rms_voltage(k) = sqrt(1/length(data(k).h_data)*sum((data(k).h_data-voltage_mean(k)).^2));
    
    Tu(k) = rms_velocity(k)/velocity_mean(k)*100;
end

figure
yyaxis left
plot(rms_velocity)
hold on
yyaxis right
plot(rms_voltage)

figure
scatter(velocity_mean,Tu)

%% Split up data for each test sequence into 0.5s segments
close all

SAMPLEFREQ = 50000;
lowpass = 1e3;
highpass = 2e1;
d1 = designfilt('lowpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',lowpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');
d2 = designfilt('highpassiir','FilterOrder',5, ...
    'HalfPowerFrequency',highpass,'SampleRate',SAMPLEFREQ,'DesignMethod','butter');


time_s = 0.5;
seglen = SAMPLEFREQ*time_s;

length_test = zeros(1,length(data));
nseg = zeros(1,length(data));
for k = 1:5
    length_test(k) = floor(length(data(k).h_v)/seglen)*seglen;
    nseg(k) = length_test(k)/seglen;
end

for k = 1:5
    data(k).mtx = zeros(nseg(k),seglen);
    data(k).seg_mean = zeros(nseg(k),1);
    
    for i = 1:nseg(k)
        data(k).mtx(i,:) = data(k).h_v((i-1)*seglen+1:i*seglen);
        data(k).seg_mean(i) = mean(data(k).mtx(i,:));
    end
end

for k = 1:5
    for i = 1:nseg(k)
        data(k).mtx(i,:) = filtfilt(d1,data(k).mtx(i,:)-data(k).seg_mean(i)) + data(k).seg_mean(i);
        data(k).mtx(i,:) = filtfilt(d2,data(k).mtx(i,:)-data(k).seg_mean(i)) + data(k).seg_mean(i);
    end
end

for k = 1:5
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
    Tu_mean(k) = mean(data(k).Tu);
    Velocity_mean(k) = mean(data(k).seg_mean);
end


figure
for k = 1:5
    scatter(data(k).seg_mean,data(k).Tu)
    hold on
end

width = 5;
height = 4.5;

scatter(Velocity_mean,Tu_mean,70,'filled','k')
ylim([0 0.7])
ylabel('Tu, $\%$','Interpreter','latex')
xlabel('Airspeed, m/s','Interpreter','latex')
set(gca,'FontSize',18)
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size









