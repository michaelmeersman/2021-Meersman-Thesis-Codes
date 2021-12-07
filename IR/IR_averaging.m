%%
clear;
close all;
clc;

width = 4.5;     % Width in inches 8
height = 5;    % Height in inches 4
alw = 1;       % AxesLineWidth 
fsz = 20;      % Fontsize 
lw = 2;        % LineWidth 
fsl = 24;      % Font size Label 
msz = 8;       % MarkerSize

%%
lims = [7500 7900]; %intensity counts range for plotting
%load mapped data
testcase = 5;
alpha = [-0.2 -0.5 2.1 2.2 3.3];
RE = [461000 513000 395000 387000 378000];
images = 1:10;

for i=1:length(images)
    setpoint = ['Mapped_Case_' num2str(testcase) '_IM' num2str(images(i)) '.mat'];
    Mapped.(['IM_' num2str(images(i))]) = load(setpoint);
end

for i=1:length(images)
    xlength(i) = length(Mapped.(['IM_' num2str(images(i))]).xc);
    zlength(i) = length(Mapped.(['IM_' num2str(images(i))]).zc);
end

[xlength,indx] = min(xlength);
[zlength,indz] = min(zlength);

xc = Mapped.(['IM_' num2str(images(indx))]).xc;
zc = Mapped.(['IM_' num2str(images(indz))]).zc;

[Z,X] = meshgrid(zc,xc);

for i=1:length(images)
    [Zaux,Xaux] = meshgrid(Mapped.(['IM_' num2str(images(i))]).zc,Mapped.(['IM_' num2str(images(i))]).xc);
    Interpolated(:,:,i) = interp2(Zaux,Xaux,double(Mapped.(['IM_' num2str(images(i))]).Area),Z,X);
end

% Time-averaging
Time = mean(Interpolated,3);
% Time-average for the non-spanwise-corrected case
Time_NC = Time;
% Time-average for the non-spanwise-corrected and without no-flow offset
Time_NC_NNF = Time;
% Time-average for no flow offset (not averaged in span)
Time_NF = Time;

%% Spanwise-correction

% Spanwise-correction using 3rd order polynomial fit to spanwise
% distribution of intensity at every x-position
for i=1:length(xc)
    third = polyfit(zc,Time(i,:),3);
    fitted(i,:) = polyval(third,zc);
end
fitted = fitted-max(fitted,[],2);
Time = Time-fitted;

%% No-flow offset

% No-flow correction to account for variation in x-direction based on an
% image taken with no flow
No_flow = load('Mapped_no_flow_50.mat');
No_flow.Interpolated = interp2(No_flow.zc,No_flow.xc,double(No_flow.Area),Z,X);
No_flow.no_SC = No_flow.Interpolated;

% images_NF = 20;
% 
% for i=1:length(images_NF)
%     setpoint = ['Mapped_no_flow_' num2str(images(i)) '.mat'];
%     No_flow.(['IM_' num2str(images(i))]) = load(setpoint);
% end
% 
% for i=1:length(images_NF)
%     [Zaux,Xaux] = meshgrid(No_flow.(['IM_' num2str(images(i))]).zc,No_flow.(['IM_' num2str(images(i))]).xc);
%     Interpolated(:,:,i) = interp2(Zaux,Xaux,double(No_flow.(['IM_' num2str(images(i))]).Area),Z,X);
% end
% 
% No_flow.Interpolated = mean(Interpolated,3);
% No_flow.no_SC = No_flow.Interpolated;
% 
% % Spanwise correction of no-flow image so that averaging makes sense
for i=1:length(xc)
    third = polyfit(zc,No_flow.Interpolated(i,:),3);
    fitted(i,:) = polyval(third,zc);
end
fitted = fitted-max(fitted,[],2);
No_flow.Interpolated = No_flow.Interpolated-fitted;

% 
% %Spanwise average of the no-flow image
No_flow_avg = mean(No_flow.Interpolated,2);

%% Normalization

% Apply no-flow offset to the spanwise-corrected map
Time_C =  -No_flow_avg + Time;
% Normalize the spanwise-corrected with no-flow offset case
cold = min(min(Time_C));
Time_C = Time_C-cold;
top = max(max(Time_C));
Time_C = Time_C/top;

% Apply no-flow offset to the non-spanwise-corrected map
Time_NC= -No_flow_avg + Time_NC;
% Normalize the non-spanwise-correction with no-flow offset case
cold = min(min(Time_NC));
Time_NC = Time_NC-cold;
top = max(max(Time_NC));
Time_NC = Time_NC/top;

% Normalize the non-spanwise-correction and without no-flow offset case
cold = min(min(Time_NC_NNF));
Time_NC_NNF = Time_NC_NNF-cold;
top = max(max(Time_NC_NNF));
Time_NC_NNF = Time_NC_NNF/top;

% Normalize the non-averaged no-flow offset case
Time_NF = Time_NF - No_flow.no_SC;
cold = min(min(Time_NF));
Time_NF = Time_NF-cold;
top = max(max(Time_NF));
Time_NF = Time_NF/top;

%% Locating separation, transition, and reattachment

% Number of points to sample for moving averages
mov_num = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No-flow offset without spanwise correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find distribution of separation, transition, and reattachment fronts
x = X(:,1);
z = Z(1,:);
reattachment_NC = zeros(size(Time_NC,2),1);
transition_NC = zeros(size(Time_NC,2),1);
separation_NC = zeros(size(Time_NC,2),1);
for i=1:size(Time_NC,2)
    T_NC = Time_NC(:,i);
    Gradient_NC = gradient(movmean(T_NC,mov_num))/(x(2)-x(1));
    %find LSB locations
    [mm,mloc] = min(T_NC);
    reattachment_NC(i) = x(mloc);
    [mg,mgloc] = min(Gradient_NC);
    transition_NC(i) = x(mgloc);
    [maxi,maxiloc] = max(Gradient_NC);
    separation_NC(i) = x(maxiloc);
end

% Find LSB locations averaged in span
T_NC = mean(Time_NC,2);
Gradient_NC = gradient(movmean(T_NC,mov_num))/(x(2)-x(1));
[mm,mloc] = min(T_NC);
rea_NC = x(mloc);
[mg,mgloc] = min(Gradient_NC);
tr_NC = x(mgloc);
[maxi,maxiloc] = max(Gradient_NC);
sep_NC = x(maxiloc);

% Filter location distribution to produce more smooth front
% separation_NC = medfilt1(separation_NC,3);
% transition_NC = medfilt1(transition_NC,3);
% reattachment_NC = medfilt1(reattachment_NC,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No-flow offset with spanwise correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find distribution of separation, transition, and reattachment fronts
x = X(:,1);
z = Z(1,:);
reattachment_C = zeros(size(Time_C,2),1);
transition_C = zeros(size(Time_C,2),1);
separation_C = zeros(size(Time_C,2),1);
for i=1:size(Time_C,2)
    T_C = Time_C(:,i);
    Gradient_C = gradient(movmean(T_C,mov_num))/(x(2)-x(1));
    %find LSB locations
    [mm,mloc] = min(T_C);
    reattachment_C(i) = x(mloc);
    [mg,mgloc] = min(Gradient_C);
    transition_C(i) = x(mgloc);
    [maxi,maxiloc] = max(Gradient_C);
    separation_C(i) = x(maxiloc);
end

% Find LSB locations averaged in span
T_C = mean(Time_C,2);
Gradient_C = gradient(movmean((T_C),mov_num))/(x(2)-x(1));
[mm,mloc] = min(T_C);
rea_C = x(mloc);
[mg,mgloc] = min(Gradient_C);
tr_C = x(mgloc);
[maxi,maxiloc] = max(Gradient_C);
sep_C = x(maxiloc);

% Filter location distribution to produce more smooth front
% separation_C = medfilt1(separation_C,3);
% transition_C = medfilt1(transition_C,3);
% reattachment_C = medfilt1(reattachment_C,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without no-flow offset and without spanwise correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find distribution of separation, transition, and reattachment fronts
x = X(:,1);
z = Z(1,:);
reattachment_NC_NNF = zeros(size(Time_NC_NNF,2),1);
transition_NC_NNF = zeros(size(Time_NC_NNF,2),1);
separation_NC_NNF = zeros(size(Time_NC_NNF,2),1);
for i=1:size(Time_NC_NNF,2)
    T_NC_NNF = Time_NC_NNF(:,i);
    Gradient_NC_NNF = gradient(movmean(T_NC_NNF,mov_num))/(x(2)-x(1));
    %find LSB locations
    [mm,mloc] = min(T_NC_NNF);
    reattachment_NC_NNF(i) = x(mloc);
    [mg,mgloc] = min(Gradient_NC_NNF);
    transition_NC_NNF(i) = x(mgloc);
    [maxi,maxiloc] = max(Gradient_NC_NNF);
    separation_NC_NNF(i) = x(maxiloc);
end

% Find LSB locations averaged in span
T_NC_NNF = mean(Time_NC_NNF,2);
Gradient_NC_NNF = gradient(T_NC_NNF)/(x(2)-x(1));
[mm,mloc] = min(T_NC_NNF);
rea_NC_NNF = x(mloc);
[mg,mgloc] = min(Gradient_NC_NNF);
tr_NC_NNF = x(mgloc);
[maxi,maxiloc] = max(Gradient_NC_NNF);
sep_NC_NNF = x(maxiloc);

% Filter location distribution to produce more smooth front
% separation_NC_NNF = medfilt1(separation_NC_NNF,3);
% transition_NC_NNF = medfilt1(transition_NC_NNF,3);
% reattachment_NC_NNF = medfilt1(reattachment_NC_NNF,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-SC no-flow offset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find distribution of separation, transition, and reattachment fronts
x = X(:,1);
z = Z(1,:);
reattachment_NF = zeros(size(Time_NF,2),1);
transition_NF = zeros(size(Time_NF,2),1);
separation_NF = zeros(size(Time_NF,2),1);
for i=1:size(Time_NF,2)
    T_NF = Time_NF(:,i);
    Gradient_NF = gradient(movmean(T_NF,mov_num))/(x(2)-x(1));
    %find LSB locations
    [mm,mloc] = min(T_NF);
    reattachment_NF(i) = x(mloc);
    [mg,mgloc] = min(Gradient_NF);
    transition_NF(i) = x(mgloc);
    [maxi,maxiloc] = max(Gradient_NF);
    separation_NF(i) = x(maxiloc);
end

% Find LSB locations averaged in span
T_NF = mean(Time_NF,2);
Gradient_NF = gradient(movmean(T_NF,mov_num))/(x(2)-x(1));
[mm,mloc] = min(T_NF);
rea_NF = x(mloc);
[mg,mgloc] = min(Gradient_NF);
tr_NF = x(mgloc);
[maxi,maxiloc] = max(Gradient_NF);
sep_NF = x(maxiloc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw no-flow offset image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cold = min(min(No_flow.no_SC));
No_flow.no_SC = No_flow.no_SC-cold;
top = max(max(No_flow.no_SC));
No_flow.no_SC = No_flow.no_SC/top;

T_raw_NF = mean(No_flow.no_SC,2);
Gradient_raw_NF = gradient(movmean(T_raw_NF,mov_num))/(x(2)-x(1));

%% Plot
close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot spanwise-corrected and no-flow offset intensity map and gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure

figure
subplot(1,2,1)
pcolor(Z,X,Time_C);   
shading('interp');
hold on
cbar = colorbar('eastoutside');
cbar.Label.String = 'Intensity';
cbar.FontSize = fsl;
cbar.Label.Interpreter = 'latex';
cbar.Ticks = [];
colormap(gca,jet(128))
set(gca,'clim',[0 1])
axis tight equal
set(gca,'YDir','reverse')
set(gca, 'Layer', 'top')
xlim([-0.3/2 0.3/2])
xlabel('z/c', 'fontsize',fsl,'interpreter','latex')
ylabel('x/c', 'fontsize',fsl,'interpreter','latex')
% Plot spanwise-averaged locations
yline(sep_C,':k','Separation','interpreter','latex','fontsize',14)
yline(tr_C,':k','Transition','interpreter','latex','fontsize',14)
yline(rea_C,':w','Reattachment','interpreter','latex','fontsize',14,'Labelverticalalignment','bottom')
% Plot median filtered, location fronts
plot3(zc,separation_C,ones(length(zc),1),'k')
plot3(zc,transition_C,ones(length(zc),1),'k')
plot3(zc,reattachment_C,ones(length(zc),1),'w')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);

subplot(1,2,2)
yyaxis left
plot(xc,T_C,'LineWidth',2)
hold on
ylabel('Intensity', 'fontsize',fsl,'interpreter','latex')
yyaxis right
plot(xc,Gradient_C,'LineWidth',2)
grid on
xline(sep_C,'--k','Separation','LabelOrientation','aligned','LabelVerticalAlignment','middle','fontsize',14,'interpreter','latex')
xline(tr_C,'--k','Transition','LabelOrientation','aligned','LabelVerticalAlignment','top','fontsize',14,'interpreter','latex')
xline(rea_C,'--k','Reattachment','LabelOrientation','aligned','LabelVerticalAlignment','top','fontsize',14,'interpreter','latex')
xlabel('x/c', 'fontsize',fsl,'interpreter','latex')
ylabel('Numerical gradient', 'fontsize',fsl,'interpreter','latex')
sgtitle(['Spanwise-correction \& No-flow offset, Re = ' num2str(RE(testcase)) ', $\alpha = ' num2str(alpha(testcase)) '^{\circ}$'],'interpreter','latex','fontsize',fsl) 
% title('Spanwise-correction \& No-flow offset', 'fontsize',fsl,'interpreter','latex')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*260, height*120]); %<- Set size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot non-spanwise-corrected and no-flow offset intensity map and gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure
figure
subplot(1,2,1)
pcolor(Z,X,Time_NC);   
shading('interp');
hold on
cbar = colorbar('eastoutside');
cbar.Label.String = 'Intensity';
cbar.FontSize = fsl;
cbar.Label.Interpreter = 'latex';
cbar.Ticks = [];
colormap(gca,jet(128))
set(gca,'clim',[0 1])
axis tight equal
set(gca,'YDir','reverse')
set(gca, 'Layer', 'top')
xlim([-0.3/2 0.3/2])
xlabel('z/c', 'fontsize',fsl,'interpreter','latex')
ylabel('x/c', 'fontsize',fsl,'interpreter','latex')
% Plot spanwise-averaged locations
yline(sep_NC,':k','Separation','interpreter','latex','fontsize',14)
yline(tr_NC,':k','Transition','interpreter','latex','fontsize',14)
yline(rea_NC,':w','Reattachment','interpreter','latex','fontsize',14)
% Plot median filtered, location fronts
plot3(zc,separation_NC,ones(length(zc),1),'k')
plot3(zc,transition_NC,ones(length(zc),1),'k')
plot3(zc,reattachment_NC,ones(length(zc),1),'w')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);

subplot(1,2,2)
yyaxis left
plot(xc,T_NC,'LineWidth',2)
hold on
ylabel('Intensity', 'fontsize',fsl,'interpreter','latex')
yyaxis right
plot(xc,Gradient_NC,'LineWidth',2)
grid on
xline(sep_NC,'--k','Separation','LabelOrientation','aligned','LabelVerticalAlignment','middle','fontsize',14,'interpreter','latex')
xline(tr_NC,'--k','Transition','LabelOrientation','aligned','LabelVerticalAlignment','top','fontsize',14,'interpreter','latex')
xline(rea_NC,'--k','Reattachment','LabelOrientation','aligned','LabelVerticalAlignment','top','fontsize',14,'interpreter','latex')
xlabel('x/c', 'fontsize',fsl,'interpreter','latex')
ylabel('Numerical gradient', 'fontsize',fsl,'interpreter','latex')
sgtitle(['No spanwise-correction, with no-flow offset, Re = ' num2str(RE(testcase)) ', $\alpha = ' num2str(alpha(testcase)) '^{\circ}$'],'interpreter','latex','fontsize',fsl) 
title('Spanwise-correction \& No-flow offset', 'fontsize',fsl,'interpreter','latex')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*260, height*120]); %<- Set size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot non-spanwise-corrected without no-flow offset intensity map and gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure
figure
subplot(1,2,1)
pcolor(Z,X,Time_NC_NNF);   
shading('interp');
hold on
cbar = colorbar('eastoutside');
cbar.Label.String = 'Intensity';
cbar.FontSize = fsl;
cbar.Label.Interpreter = 'latex';
cbar.Ticks = [];
colormap(gca,jet(128))
set(gca,'clim',[0 1])
axis tight equal
set(gca,'YDir','reverse')
set(gca, 'Layer', 'top')
xlim([-0.3/2 0.3/2])
xlabel('z/c', 'fontsize',fsl,'interpreter','latex')
ylabel('x/c', 'fontsize',fsl,'interpreter','latex')
% Plot spanwise-averaged locations
yline(sep_NC_NNF,':k','Separation','interpreter','latex','fontsize',14)
yline(tr_NC_NNF,':k','Transition','interpreter','latex','fontsize',14)
yline(rea_NC_NNF,':w','Reattachment','interpreter','latex','fontsize',14)
% Plot median filtered, location fronts
plot3(zc,separation_NC_NNF,ones(length(zc),1),'k')
plot3(zc,transition_NC_NNF,ones(length(zc),1),'k')
plot3(zc,reattachment_NC_NNF,ones(length(zc),1),'w')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);

subplot(1,2,2)
yyaxis left
plot(xc,T_NC_NNF,'LineWidth',2)
hold on
ylabel('Intensity', 'fontsize',fsl,'interpreter','latex')
yyaxis right
plot(xc,Gradient_NC_NNF,'LineWidth',2)
grid on
xline(sep_NC_NNF,'--k','Separation','LabelOrientation','aligned','LabelVerticalAlignment','middle','fontsize',14,'interpreter','latex')
xline(tr_NC_NNF,'--k','Transition','LabelOrientation','aligned','LabelVerticalAlignment','top','fontsize',14,'interpreter','latex')
xline(rea_NC_NNF,'--k','Reattachment','LabelOrientation','aligned','LabelVerticalAlignment','top','fontsize',14,'interpreter','latex')
xlabel('x/c', 'fontsize',fsl,'interpreter','latex')
ylabel('Numerical gradient', 'fontsize',fsl,'interpreter','latex')
sgtitle(['No spanwise-correction, without no-flow offset, Re = ' num2str(RE(testcase)) ', $\alpha = ' num2str(alpha(testcase)) '^{\circ}$'],'interpreter','latex','fontsize',fsl) 
% title('Spanwise-correction \& No-flow offset', 'fontsize',fsl,'interpreter','latex')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*260, height*120]); %<- Set size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot non-SC no-flow offset case
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure
figure
subplot(1,2,1)
pcolor(Z,X,Time_NF);   
shading('interp');
hold on
cbar = colorbar('eastoutside');
cbar.Label.String = 'Intensity';
cbar.FontSize = fsl;
cbar.Label.Interpreter = 'latex';
cbar.Ticks = [];
colormap(gca,jet(128))
set(gca,'clim',[0 1])
axis tight equal
set(gca,'YDir','reverse')
set(gca, 'Layer', 'top')
xlim([-0.3/2 0.3/2])
xlabel('z/c', 'fontsize',fsl,'interpreter','latex')
ylabel('x/c', 'fontsize',fsl,'interpreter','latex')
% Plot spanwise-averaged locations
yline(sep_NF,':k','Separation','interpreter','latex','fontsize',14)
yline(tr_NF,':k','Transition','interpreter','latex','fontsize',14)
yline(rea_NF,':w','Reattachment','interpreter','latex','fontsize',14)
% Plot median filtered, location fronts
plot3(zc,separation_NF,ones(length(zc),1),'k')
plot3(zc,transition_NF,ones(length(zc),1),'k')
plot3(zc,reattachment_NF,ones(length(zc),1),'w')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);

subplot(1,2,2)
yyaxis left
plot(xc,T_NF,'LineWidth',2)
hold on
ylabel('Intensity', 'fontsize',fsl,'interpreter','latex')
yyaxis right
plot(xc,Gradient_NF,'LineWidth',2)
grid on
xline(sep_NF,'--k','Separation','LabelOrientation','aligned','LabelVerticalAlignment','middle','fontsize',14,'interpreter','latex')
xline(tr_NF,'--k','Transition','LabelOrientation','aligned','LabelVerticalAlignment','top','fontsize',14,'interpreter','latex')
xline(rea_NF,'--k','Reattachment','LabelOrientation','aligned','LabelVerticalAlignment','top','fontsize',14,'interpreter','latex')
xlabel('x/c', 'fontsize',fsl,'interpreter','latex')
ylabel('Numerical gradient', 'fontsize',fsl,'interpreter','latex')
sgtitle(['Raw no-flow offset, Re = ' num2str(RE(testcase)) ', $\alpha = ' num2str(alpha(testcase)) '^{\circ}$'],'interpreter','latex','fontsize',fsl) 
title('Spanwise-correction \& No-flow offset', 'fontsize',fsl,'interpreter','latex')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*260, height*120]); %<- Set size


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot just the raw no-flow offset image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%figure
figure
subplot(1,2,1)
pcolor(zc,xc,No_flow.no_SC);   
shading('interp');
hold on
cbar = colorbar('eastoutside');
cbar.Label.String = 'Intensity';
cbar.FontSize = fsl;
cbar.Label.Interpreter = 'latex';
cbar.Ticks = [];
colormap(gca,jet(128))
set(gca,'clim',[0 1])
axis tight equal
set(gca,'YDir','reverse')
set(gca, 'Layer', 'top')
xlim([-0.3/2 0.3/2])
xlabel('z/c', 'fontsize',fsl,'interpreter','latex')
ylabel('x/c', 'fontsize',fsl,'interpreter','latex')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);

subplot(1,2,2)
yyaxis left
plot(xc,T_raw_NF,'LineWidth',2)
hold on
ylabel('Intensity', 'fontsize',fsl,'interpreter','latex')
yyaxis right
plot(xc,Gradient_raw_NF,'LineWidth',2)
grid on
xlabel('x/c', 'fontsize',fsl,'interpreter','latex')
ylabel('Numerical gradient', 'fontsize',fsl,'interpreter','latex')
sgtitle('Raw no-flow Image','interpreter','latex','fontsize',fsl) 
% title('Spanwise-correction \& No-flow offset', 'fontsize',fsl,'interpreter','latex')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*260, height*120]); %<- Set size
% save(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\corrected_case_' num2str(testcase) '.mat'])