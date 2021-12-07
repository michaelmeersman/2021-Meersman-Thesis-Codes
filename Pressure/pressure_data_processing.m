clear; close all

% Import pressure tap locations from file
taps = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\' ...
     '7. Matlab Pressure Tap Codes\6series_pressure_tap_FLIGHT.csv']);
upper_mid_taps = taps(1:30,1);
lower_mid_taps = taps(1:8,2);

% Flight Model Pressure tap indices
up_in = 1:4;
up_mid = 5:34;
up_out = 35:38;
low_in = 39:41;
low_mid = 42:49;
low_out = 50:52;

%%% Wind Tunnel Model Pressure tap indices
top1 = 1:7;
top2 = 8:34;
top3 = 35:41;
bottom1 = 42:46;
bottom2 = 47:55;
bottom3 = 56:60;
static = 61;
total = 62;

% Load Calibration/sensor offset file
flight_no = 3;
cal = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '1.1 Flight Test 9_18_21\RPI_Flight_' num2str(flight_no) '\CALIBRATION_09_18_21_Flight_' num2str(flight_no) '.txt'],'\t');

% Test case consistent regions
eval_sect = [1 497;
             1 1096;
             141 962;
             457 880;
             1180 2160];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Flight Test data file%%%%%%%%%%%%%%%%%%%
testcase = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
daq = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '1.1 Flight Test 9_18_21\RPI_Flight_' num2str(flight_no) '\pressure_data_case_' num2str(testcase) '_f.csv']);
alldata = daq(:,:);

% Parse data types
p = alldata(:,5:end);
alpha = alldata(:,1);
pitot = alldata(:,2);
gps = alldata(:,3);

% Apply calibration to pressure data and convert to inH2O
for i = 1:length(p)
    p(i,:) = (p(i,:) - cal)./1000;
end

% q = mean(abs(p(:,end)));             % Dynamic pressure, inH2O
alpha_mean = mean(alpha);       % Average AOA, degrees
v_ind = mean(pitot);            % Average indicated airspeed, m/s
v_g = mean(gps);                % Average ground speed, m/s
mu = 1.894e-5;                  % Dynamic viscosity, m/s^2
rho = 1.13;                     % Air density, kg/m^3
c = 0.305;                      % Chord length, m
v_comp = (v_ind+v_g)/2;
q = 0.5*rho*v_comp^2/248.84;
Re = rho*v_comp*cosd(35)*c/mu;   % Reynolds number based on indicated 
                                % airspeed in chordwise direction

% Calculate Cp based on chordwise dynamic pressure
% Pressure sensors are already referenced to static pressure
Cp = zeros(length(p),52);
for i = 1:length(p)
    Cp(i,:) = (p(i,:))./q/cosd(35)^2;
end

% Extract different sections of pressure taps
p_upper_mid = mean(p(:,up_mid));
cp_upper_mid = mean(Cp(:,up_mid));
p_lower_mid = mean(p(:,low_mid));
cp_lower_mid = mean(Cp(:,low_mid));

% Calculate Std for error bars
cp_std = zeros(52,1);
for i = 1:52
    cp_std(i) = std(Cp(:,i));
end

std_upper_mid = cp_std(up_mid);
std_lower_mid = cp_std(low_mid);

% Import Wind Tunnel and XFOIL data based on the flight case
angle_offset = 3; % Angle used to match exp. Cp with Xfoil
alpha = alpha - angle_offset;
alpha_xfoil = round(alpha_mean-angle_offset,1);
up = ceil(alpha_xfoil * 4) / 4;
down = floor(alpha_xfoil * 4) / 4;
if abs(alpha_xfoil - up) >= abs(alpha_xfoil - down)
    alpha_r = down;
elseif abs(alpha_xfoil - up) < abs(alpha_xfoil - down)
    alpha_r = up;
end
    
if alpha_r < 0 && round(alpha_r) ~= alpha_r
    alpha_wt = ['m' num2str(abs(ceil(alpha_r))) 'p' num2str(100*abs(alpha_r-ceil(alpha_r)))];
elseif alpha_r < 0 && round(alpha_r) == alpha_r
    alpha_wt = ['m' num2str(abs(alpha_r))];
elseif alpha_r > 0 && round(alpha_r) == alpha_r
    alpha_wt = num2str(alpha_r);
elseif alpha_r > 0 && round(alpha_r) ~= alpha_r
    alpha_wt = [num2str(floor(alpha_r)) 'p' num2str(100*abs(alpha_r-floor(alpha_r)))];
elseif alpha_r == 0
    alpha_wt = num2str(0);
end
Re_r = round(Re,-4);
xpos = readmatrix(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '2.1 All WT Data\xpos.csv']);
WT_x_c = xpos/12;
WT = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '2.1 All WT Data\Taps_AoA_'...
    alpha_wt ...
    '_Re_200k.CSV'],',');
WT_p = WT.data(:,3:62);
WT_p([49 50],:) = [];
WT_static = WT.data(:,64);
WT_total = WT.data(:,65);
WT_q = 0.5*1.11*10.8^2;
WT_Cp = (mean(WT_p)-mean(WT_static))*6894.76/WT_q;

%%%%%%%%%%% Include 300k WT Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_mean = round(alpha_mean-angle_offset) - 1;
if  a_mean < 0
        aoastr = ['n' num2str(abs(a_mean))];
elseif a_mean >= 0
        aoastr = num2str(a_mean);
end
data_300k = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA'...
        '\Summer WT Data 300k\07_06_300k_' aoastr '.txt'],'\t',4);
WT_Cp_300 = data_300k.data(:,8);
WT_x_c_300 = data_300k.data(:,3)/12;

wing = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '2.1 All WT Data\Modified_NACA _643_618.dat']);
wing = wing.data;
wing([44 87],:) = [];
[polar, solution] = xfoil(wing,alpha_xfoil,Re_r,0,'ppar n 200','oper/vpar n 10','oper/iter 500');
xfoil_x_c = solution.xcp;
xfoil_Cp = solution.cp;
x_close = [upper_mid_taps(end) lower_mid_taps(end)];
cp_close = [cp_upper_mid(end) cp_lower_mid(end)];

% Error bars for pressure sensors
err_upper = 0.0025*20/q+std_upper_mid;
err_lower = 0.0025*20/q+std_lower_mid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width = 6.5;
height = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create plot of mean pressure distribution comparison between fligh, xflr and wind tunnel
h(1) = figure;
subplot(7,1,[1:4])
plot(upper_mid_taps,cp_upper_mid,'k','LineWidth',1.5)
hold on
plot(xfoil_x_c,xfoil_Cp,'-.b','LineWidth',1.5)
plot(WT_x_c(horzcat(top2,fliplr(bottom2))),WT_Cp(horzcat(top2,fliplr(bottom2))),'--r','LineWidth',1.5)
plot(WT_x_c_300(horzcat(top2,fliplr(bottom2))),WT_Cp_300(horzcat(top2,fliplr(bottom2))),'--','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
scatter(upper_mid_taps,cp_upper_mid,30,'k','filled')
plot(lower_mid_taps,cp_lower_mid,'k','LineWidth',1.5)
scatter(lower_mid_taps,cp_lower_mid,30,'k','filled')
plot(x_close,cp_close,'k','LineWidth',1.5)
errorbar(upper_mid_taps,cp_upper_mid,err_upper,'k','LineStyle','none')
errorbar(lower_mid_taps,cp_lower_mid,err_lower,'k','LineStyle','none')
ylim([min([cp_upper_mid cp_lower_mid])-0.2 max([cp_upper_mid cp_lower_mid])+0.25])
xlim([-0.05 1])
legend(['Flight, $Re_{c} = $' num2str(round(Re,-3), '%.f') ', $\alpha = $' num2str(alpha_xfoil)],...
       ['XFOIL, $Re = $' num2str(Re_r) ', $\alpha = $' num2str(alpha_xfoil)],...
       ['ALSWT Unswept, $Re = 200k$, $\alpha = $' num2str(alpha_r)],...
       ['ALSWT Unswept, $Re = 300k$, $\alpha = $' num2str(round(alpha_r)) ],...
       'location','south','FontSize',14,'interpreter','latex')
set(gca, 'YDir','reverse')
set(gca, 'LineWidth', 1.5, 'FontSize', 15)
xlabel('$x/c$','Interpreter','latex','FontSize',24)
ylabel('$C_p$','Interpreter','latex','FontSize',24)
% set(gca, 'FontName','Times','FontWeight','bold')
ax.LabelFontSizeMultiplier = 20;
% title(['NACA $64_{3}-618$ Flight Test $C_p$ (' num2str(v_ind,2) 'm/s Indicated)'],'Interpreter','latex','FontSize',24)
title([ '$Re = $' num2str(Re_r) ', $\alpha = $' num2str(alpha_xfoil)],'Interpreter','latex','FontSize',24)
grid on
hold off


% Time-trace plot of AOA and airspeed
subplot(7,1,[6:7])
time = (0:length(alpha)-1)/130;
yyaxis left
plot(time,alpha,'Linewidth',1)
ylab1 = ylabel('$\alpha$ $^o$','interpreter','latex');
ylim([min(alpha)-1.5 max(alpha)+1.5])
hold on
yyaxis right
plot(time,pitot,'Linewidth',1)
ylab2 = ylabel('Velocity, $m/s$','interpreter','latex');
xlim([0 (length(alpha)+2)/130])
xticks([0:2:18])
ylim([min(pitot)-3 max(pitot)+3])
xlab = xlabel('Time, s','interpreter','latex');
% title('Angle of Attack \& Airspeed Time-traces','interpreter','latex')
set(gca, 'FontSize',16)
ylab1.FontSize = 24;
ylab2.FontSize = 22;
xlab.FontSize = 24;
% Show Selected range for conditional filtering
plot([eval_sect(testcase,1)/130 eval_sect(testcase,2)/130 eval_sect(testcase,2)/130 eval_sect(testcase,1)/130 eval_sect(testcase,1)/130],...
    [max(pitot)+ 2 max(pitot)+ 2 min(pitot)- 2 min(pitot)- 2 max(pitot)+ 2],'r','LineWidth',1.5)

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-500 width*100, height*150]); %<- Set size


%% Second plot with only selected data region

alldata = daq([eval_sect(testcase,1):eval_sect(testcase,2)],:);

% Parse data types
p = alldata(:,5:end);
alpha = alldata(:,1);
pitot = alldata(:,2);
gps = alldata(:,3);

% Apply calibration to pressure data and convert to inH2O
for i = 1:length(p)
    p(i,:) = (p(i,:) - cal)./1000;
end


alpha_mean = mean(alpha);       % Average AOA, degrees
v_ind = mean(pitot);            % Average indicated airspeed, m/s
% q = mean(abs(p(:,end)));             % Dynamic pressure, inH2O
v_g = mean(gps);                % Average ground speed, m/s
mu = 1.894e-5;                  % Dynamic viscosity, m/s^2
rho = 1.13;                     % Air density, kg/m^3
c = 0.305;                      % Chord length, m
v_comp = (v_ind+v_g)/2;
q = 0.5*rho*v_comp^2/248.84;
Re = rho*v_comp*cosd(35)*c/mu;   % Reynolds number based on indicated 
                                % airspeed in chordwise direction

% Calculate Cp based on chordwise dynamic pressure
% Pressure sensors are already referenced to static pressure
Cp = zeros(length(p),52);
for i = 1:length(p)
    Cp(i,:) = (p(i,:))./q/cosd(35)^2;
end

% Extract different sections of pressure taps
p_upper_mid = mean(p(:,up_mid));
cp_upper_mid = mean(Cp(:,up_mid));
p_lower_mid = mean(p(:,low_mid));
cp_lower_mid = mean(Cp(:,low_mid));

% Calculate Std for error bars
cp_std = zeros(52,1);
for i = 1:52
    cp_std(i) = std(Cp(:,i));
end

std_upper_mid = cp_std(up_mid);
std_lower_mid = cp_std(low_mid);

% Import Wind Tunnel and XFOIL data based on the flight case
angle_offset = 3; % Angle used to match exp. Cp with Xfoil
alpha = alpha - angle_offset;
alpha_xfoil = round(alpha_mean-angle_offset,1);
up = ceil(alpha_xfoil * 4) / 4;
down = floor(alpha_xfoil * 4) / 4;
if abs(alpha_xfoil - up) >= abs(alpha_xfoil - down)
    alpha_r = down;
elseif abs(alpha_xfoil - up) < abs(alpha_xfoil - down)
    alpha_r = up;
end
    
if alpha_r < 0 && round(alpha_r) ~= alpha_r
    alpha_wt = ['m' num2str(abs(ceil(alpha_r))) 'p' num2str(100*abs(alpha_r-ceil(alpha_r)))];
elseif alpha_r < 0 && round(alpha_r) == alpha_r
    alpha_wt = ['m' num2str(abs(alpha_r))];
elseif alpha_r > 0 && round(alpha_r) == alpha_r
    alpha_wt = num2str(alpha_r);
elseif alpha_r > 0 && round(alpha_r) ~= alpha_r
    alpha_wt = [num2str(floor(alpha_r)) 'p' num2str(100*abs(alpha_r-floor(alpha_r)))];
elseif alpha_r == 0
    alpha_wt = num2str(0);
end
Re_r = round(Re,-3);
xpos = readmatrix(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '2.1 All WT Data\xpos.csv']);
WT_x_c = xpos/12;
WT = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '2.1 All WT Data\Taps_AoA_'...
    alpha_wt ...
    '_Re_200k.CSV'],',');
WT_p = WT.data(:,3:62);
WT_p([49 50],:) = [];
WT_static = WT.data(:,64);
WT_total = WT.data(:,65);
WT_q = 0.5*1.11*10.8^2;
WT_Cp = (mean(WT_p)-mean(WT_static))*6894.76/WT_q;

%%%%%%%%%%% Include 300k WT Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_mean = round(alpha_mean-angle_offset) - 1;
if  a_mean < 0
        aoastr = ['n' num2str(abs(a_mean))];
elseif a_mean >= 0
        aoastr = num2str(a_mean);
end
data_300k = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA'...
        '\Summer WT Data 300k\07_06_300k_' aoastr '.txt'],'\t',4);
WT_Cp_300 = data_300k.data(:,8);
WT_x_c_300 = data_300k.data(:,3)/12;

%%%%%%%%%%%%% XFOIL Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wing = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '2.1 All WT Data\Modified_NACA _643_618.dat']);
wing = wing.data;
wing([44 87],:) = [];
[polar, solution] = xfoil(wing,alpha_xfoil,Re_r,0,'ppar n 200','oper/vpar n 7','oper/iter 500');
xfoil_x_c = solution.xcp;
xfoil_Cp = solution.cp;
x_close = [upper_mid_taps(end) lower_mid_taps(end)];
cp_close = [cp_upper_mid(end) cp_lower_mid(end)];

% Error bars for pressure sensors
err_upper = 0.0025*20/q+std_upper_mid;
err_lower = 0.0025*20/q+std_lower_mid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width = 6.5;
height = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create plot of mean pressure distribution comparison between fligh, xflr and wind tunnel
h(2) = figure;
subplot(7,1,[1:4])
plot(upper_mid_taps,cp_upper_mid,'k','LineWidth',1.5)
hold on
plot(xfoil_x_c,xfoil_Cp,'-.b','LineWidth',1.5)
plot(WT_x_c(horzcat(top2,fliplr(bottom2))),WT_Cp(horzcat(top2,fliplr(bottom2))),'--r','LineWidth',1.5)
plot(WT_x_c_300(horzcat(top2,fliplr(bottom2))),WT_Cp_300(horzcat(top2,fliplr(bottom2))),'--','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
scatter(upper_mid_taps,cp_upper_mid,30,'k','filled')
plot(lower_mid_taps,cp_lower_mid,'k','LineWidth',1.5)
scatter(lower_mid_taps,cp_lower_mid,30,'k','filled')
plot(x_close,cp_close,'k','LineWidth',1.5)
errorbar(upper_mid_taps,cp_upper_mid,err_upper,'k','LineStyle','none')
errorbar(lower_mid_taps,cp_lower_mid,err_lower,'k','LineStyle','none')
ylim([min([cp_upper_mid cp_lower_mid])-0.2 max([cp_upper_mid cp_lower_mid])+0.25])
xlim([-0.05 1])
legend(['Flight, $Re_{c} = $' num2str(round(Re,-3), '%.f') ', $\alpha = $' num2str(alpha_xfoil)],...
       ['XFOIL, $Re = $' num2str(Re_r) ', $\alpha = $' num2str(alpha_xfoil)],...
       ['ALSWT Unswept, $Re = 200k$, $\alpha = $' num2str(alpha_r)],...
       ['ALSWT Unswept, $Re = 300k$, $\alpha = $' num2str(round(alpha_r)) ],...
       'location','south','FontSize',14,'interpreter','latex')
set(gca, 'YDir','reverse')
set(gca, 'LineWidth', 1.5, 'FontSize', 15)
xlabel('$x/c$','Interpreter','latex','FontSize',24)
ylabel('$C_p$','Interpreter','latex','FontSize',24)
% set(gca, 'FontName','Times','FontWeight','bold')
% ax.LabelFontSizeMultiplier = 20;
% title(['NACA $64_{3}-618$ Flight Test $C_p$ (' num2str(v_ind,2) 'm/s Indicated)'],'Interpreter','latex','FontSize',24)
title([ '$Re = $' num2str(Re_r) ', $\alpha = $' num2str(alpha_xfoil)],'Interpreter','latex','FontSize',24)
grid on
hold off


% Time-trace plot of AOA and airspeed
subplot(7,1,[6:7])
time = (0:length(alpha)-1)/130;
yyaxis left
plot(time,alpha,'Linewidth',1)
ylab1 = ylabel('$\alpha$ $^o$','interpreter','latex');
ylim([min(alpha)-1.5 max(alpha)+1.5])
hold on
yyaxis right
plot(time,pitot,'Linewidth',1)
ylab2 = ylabel('Velocity, $m/s$','interpreter','latex');
xlim([-1/130 (length(alpha)+2)/130])
xticks([0:1:7])
ylim([min(pitot)-3 max(pitot)+3])
xlab = xlabel('Time, s','interpreter','latex');
% title('Angle of Attack \& Airspeed Time-traces','interpreter','latex')
set(gca, 'FontSize',16)
ylab1.FontSize = 24;
ylab2.FontSize = 22;
xlab.FontSize = 24;
% Show Selected range for conditional filtering
plot([0 (length(alpha)-1)/130 (length(alpha)-1)/130 0 0],...
    [max(pitot)+ 2 max(pitot)+ 2 min(pitot)- 2 min(pitot)- 2 max(pitot)+ 2],'r','LineWidth',1.5)

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-500 width*100, height*150]); %<- Set size
%%
% Plot including 300k wind tunnel case
h(3) = figure;
plot(upper_mid_taps,cp_upper_mid,'k','LineWidth',1.5)
hold on
plot(xfoil_x_c,xfoil_Cp,'-.b','LineWidth',1.5)
plot(WT_x_c(horzcat(top2,fliplr(bottom2))),WT_Cp(horzcat(top2,fliplr(bottom2))),'--r','LineWidth',1.5)
plot(WT_x_c_300(horzcat(top2,fliplr(bottom2))),WT_Cp_300(horzcat(top2,fliplr(bottom2))),'--','LineWidth',1.5,'Color',[0.4940 0.1840 0.5560])
scatter(upper_mid_taps,cp_upper_mid,30,'k','filled')
plot(lower_mid_taps,cp_lower_mid,'k','LineWidth',1.5)
scatter(lower_mid_taps,cp_lower_mid,30,'k','filled')
plot(x_close,cp_close,'k','LineWidth',1.5)
errorbar(upper_mid_taps,cp_upper_mid,err_upper,'k','LineStyle','none')
errorbar(lower_mid_taps,cp_lower_mid,err_lower,'k','LineStyle','none')
ylim([min([cp_upper_mid cp_lower_mid])-0.2 max([cp_upper_mid cp_lower_mid])])
xlim([-0.05 1])
legend(['Flight, $Re_{c} = $' num2str(round(Re,-3), '%.f') ', $\alpha =$ ' num2str(alpha_xfoil)],...
       ['XFOIL, $Re = $' num2str(Re_r) ', $\alpha = $' num2str(alpha_xfoil)],...
       ['ALSWT Unswept, $Re = 200k$, $\alpha = $' num2str(alpha_r)],...
       ['ALSWT Unswept, $Re = 300k$, $\alpha = $' num2str(round(alpha_r)) ],...
       'Interpreter','latex', 'FontSize',16,'location','south')
set(gca, 'YDir','reverse')
set(gca, 'LineWidth', 1.5, 'FontSize', 15)
xlab = xlabel('$x/c$','Interpreter','latex');
ylab = ylabel('$C_p$','Interpreter','latex');
set(gca, 'Fontsize',16)
xlab.FontSize = 24;
ylab.FontSize = 24;
title([ '$Re = $' num2str(Re_r) ', $\alpha = $' num2str(alpha_xfoil)],'Interpreter','latex','FontSize',24)
grid on
width = 6.3;
height =5;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-500 width*120, height*120]); %<- Set size
hold off

save(['cp_selected_case_' num2str(testcase) '.mat'])
%% Overlay IR image with Cp
% close all
load(['corrected_case_' num2str(testcase) '.mat']);

width = 2.8;
height = 4.2;

N = 1;
Gradient = [Gradient_C Gradient_NC Gradient_NC_NNF];
sep = [sep_C sep_NC sep_NC_NNF];
tr = [tr_C tr_NC tr_NC_NNF];
rea = [rea_C rea_NC rea_NC_NNF];

h(4) = figure;
% Optionally include colormap alongside the Cp plot
% subplot(1,8,1:3)
% pcolor(Z,X,Time_C);   
% shading('interp');
% hold on
% cbar = colorbar('eastoutside');
% cbar.Label.String = 'Intensity';
% cbar.FontSize = fsl;
% cbar.Label.Interpreter = 'latex';
% cbar.Ticks = [];
% colormap(gca,jet(128))
% set(gca,'clim',[0 1])
% axis tight
% set(gca,'YDir','reverse')
% set(gca, 'Layer', 'top')
% xlim([-0.3/2 0.3/2])
% % title('Spanwise-correction \& No-flow offset', 'fontsize',fsl,'interpreter','latex')
% xlabel('z/c', 'fontsize',fsl,'interpreter','latex')
% ylabel('x/c', 'fontsize',fsl,'interpreter','latex')
% set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
% % yline(mean(SEP),'k','Separation','interpreter','latex','FontSize',14)
% % yline(mean(TRANS),'k','Transition','interpreter','latex','FontSize',14)
% % yline(mean(REATT),'w','Reattachment','interpreter','latex','FontSize',14)
% yline(sep_C,':k','Separation','interpreter','latex','fontsize',14)
% yline(tr_C,':k','Transition','interpreter','latex','fontsize',14)
% yline(rea_C,':w','Reattachment','interpreter','latex','fontsize',14,'Labelverticalalignment','bottom')
% plot3(zc,separation_C,ones(length(zc),1),'k')
% plot3(zc,transition_C,ones(length(zc),1),'k')
% plot3(zc,reattachment_C,ones(length(zc),1),'w')
% set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);

% subplot(1,8,5:8)
yyaxis left
plot(upper_mid_taps,cp_upper_mid,'k','LineWidth',1.5)
hold on
upm = scatter(upper_mid_taps,cp_upper_mid,30,'k','filled');
upm.Annotation.LegendInformation.IconDisplayStyle = 'off';
lwm = plot(lower_mid_taps,cp_lower_mid,'-k','LineWidth',1.5);
lwm.Annotation.LegendInformation.IconDisplayStyle = 'off';
lwmdc = scatter(lower_mid_taps,cp_lower_mid,30,'k','filled');
lwmdc.Annotation.LegendInformation.IconDisplayStyle = 'off';
clse = plot(x_close,cp_close,'k','LineWidth',1.5,'Marker','none');
clse.Annotation.LegendInformation.IconDisplayStyle = 'off';
set(gca,'YDir','reverse')
ylab1 = ylabel('$C_p$','interpreter','latex');
ylim([min(cp_upper_mid)-0.3 max([cp_lower_mid cp_upper_mid])])
ax=gca;
ax.YAxis(1).Color = 'k';


yyaxis right
plot(xc,Gradient(:,N),'LineWidth',1.5)
grid on

%%%%%% Lines and text to indication key locations
xline(sep(N),'--k','LineWidth',1.5)
text(sep(N),1.1,'Sep','horizontalalignment','right','VerticalAlignment','bottom','Rotation',90,'interpreter','latex','FontSize',16)
xline(tr(N),'--k','LineWidth',1.5)
text(tr(N),0.5,'Trans','horizontalalignment','right','VerticalAlignment','top','Rotation',90,'interpreter','latex','FontSize',16)
xline(rea(N),'--k','LineWidth',1.5)
text(rea(N),-2,'Reatt','horizontalalignment','right','VerticalAlignment','top','Rotation',90,'interpreter','latex','FontSize',16)
%%%%%% Create box indicating region of IR camera
irbox = [0.4    max(Gradient(:,N))+1;
         0.84   max(Gradient(:,N))+1;
         0.84   min(Gradient(:,N))-1;
         0.4    min(Gradient(:,N))-1;
         0.4    max(Gradient(:,N))+1];
plot(irbox(:,1),irbox(:,2),':','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5)
text(0.85,-2,'Camera Viewing Region','Color',[0.6350 0.0780 0.1840],'LineWidth',10,'interpreter','latex','Rotation',90,'VerticalAlignment','top','FontSize',14)
xlab2 = xlabel('x/c','interpreter','latex');
ylabel('Gradient', 'fontsize',26,'interpreter','latex')
ylim([min(Gradient(:,N))-2 max(Gradient(:,N))+2])
title([ '$Re = $' num2str(Re_r) ', $\alpha = $' num2str(alpha_xfoil)],'Interpreter','latex','FontSize',24)
legend(['Flight $C_p$'],...
       ['IR Gradient'],'FontSize',14,'interpreter','latex','location', 'northwest')
% sgtitle(['Re = ' num2str(Re_r) ', $\alpha = ' num2str(alpha_xfoil) '^{\circ}$'],'interpreter','latex','fontsize',fsl) 
% title('Spanwise-correction \& No-flow offset', 'fontsize',fsl,'interpreter','latex')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
box on
ylab1.FontSize = fsz +5;
xlab2.FontSize = fsz +5;
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*260, height*150]); %<- Set size


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE SAVING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\Figures\'...
    'Cp\'];
fig1 = ['cp_time_traces_unfiltered_' num2str(testcase)];
fig2 = ['cp_time_traces_filtered_' num2str(testcase)];
fig3 = ['cp_filtered_' num2str(testcase)];
fig4 = ['IR_cp_' num2str(testcase)];

exportgraphics(h(1),[filepath fig1 '.png'],'Resolution',1000);
exportgraphics(h(2),[filepath fig2 '.png'],'Resolution',1000);
exportgraphics(h(3),[filepath fig3 '.png'],'Resolution',1000);
exportgraphics(h(4),[filepath fig4 '.png'],'Resolution',1000);

exportgraphics(h(1),[filepath fig1 '.eps']);
exportgraphics(h(2),[filepath fig2 '.eps']);
exportgraphics(h(3),[filepath fig3 '.eps']);
exportgraphics(h(4),[filepath fig4 '.eps']);

savefig(h(1),[filepath fig1 '.fig']);
savefig(h(2),[filepath fig2 '.fig']);
savefig(h(3),[filepath fig3 '.fig']);
savefig(h(4),[filepath fig4 '.fig']);






