%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  IR manual image mapping %%%%%%%%%%%%%%%%%%%%%%%%
%%%% Adrian Grille Guerra - UArizona - September 2021 %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

width = 6;     % Width in inches 8
height = 5;    % Height in inches 4
alw = 1;       % AxesLineWidth 
fsz = 20;      % Fontsize 
lw = 2;        % LineWidth 
fsl = 24;      % Font size Label 
msz = 8;       % MarkerSize
%%
lims = [7900 8500]; %intensity counts range for plotting
ROI = [151 400; 252 500]; %contains the region of interest
Im_number = 50;
Im = imread('20210918_105310_SEQ_A.tiff',Im_number);
Image = Im(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));

%marker locations
%x/c = 0.40, 0.50, 0.66, 0.84 (chord is 12 inches)
%spanwise distance 4.6 inches (roughly 100 pixels)
%markers in pixels: 73, 99, 141, 188

%create dummy of fixed image with control points
fP = [75 73; 175 73; 75 99; 175 99; 75 141; 175 141; 75 188; 175 188];
fixed = ones(size(Image));
for i=1:length(fP(:,1))
    fixed(fP(i,2)-5:fP(i,2)+5,fP(i,1)-5:fP(i,1)+5) = 0;
end

%normalize image for control point selection tool
moving = Image-min(min(Image));
moving = double(moving);
moving = moving/max(max(moving));
[MovingPoints,FixedPoints] = cpselect(moving,fixed,fP,fP,'Wait',true);
%check selection is valid
if isequal(MovingPoints,FixedPoints)
    return %there was not control point selection so it stops here
end

%projective transformation
tform = fitgeotrans(MovingPoints,FixedPoints,'projective');
Dewarped = imwarp(Image,tform);

%normalize dewarped image for control point selection tool
moving = Dewarped-min(min(Image));
moving(moving<0) = 0;
moving = double(moving);
moving = moving/max(max(moving));
%offset based on preliminary image
fPP(:,1) = fP(:,1)+10;
fPP(:,2) = fP(:,2)+40;
[MovingPoints,FixedPoints] = cpselect(moving,fixed,fPP,fP,'Wait',true);

%check which markers were detected
Markers = find(ismember(fP,FixedPoints,'rows')); 
%left location from odd markers
odd = find(mod(Markers,2));
left = mean(MovingPoints(odd,1));
%right location from even markers
even = find(~mod(Markers,2));
right = mean(MovingPoints(even,1));
%upstream location from earliest marker
earliest = Markers(1);
upstream = MovingPoints(1,2);
switch earliest
    case {1,2}
        x_upstream = 0.4;
    case {3,4}
        x_upstream = 0.5;
    otherwise
        return
end
%downstream location from latest marker
latest = Markers(end);
downstream = MovingPoints(end,2);
switch latest
    case {7,8}
        x_downstream = 0.84;
    case {5,6}
        x_downstream = 0.66;
    otherwise
        return
end
%retain only 3.6 inches (out of 4.6), 30% chord, of spanwise distance to avoid markers
tot_distance = right-left;
distance = tot_distance*3.6/4.6;
gap = (tot_distance-distance)/2;
left = round(left+gap);
right = round(right-gap);
%extract region
region = [left upstream right-left downstream-upstream];
Area = imcrop(Dewarped,region);
xc = linspace(x_upstream,x_downstream,length(Area(:,1))); %chordwise location
zc = linspace(-0.3/2,0.3/2,length(Area(1,:))); %spanwise location

%% Plots
[Z,X] = meshgrid(zc,xc); %grid for plot

figure
% subplot(1,2,1)
imagesc(Dewarped)
colormap(gca,gray)
cbar = colorbar;
cbar.Label.String = 'Intensity (counts)';
cbar.FontSize = fsz;
cbar.FontName = 'Times-Roman';
cbar.Location = 'eastoutside';
set(gca,'clim',lims)
axis tight equal
hold on
plot(MovingPoints(:,1),MovingPoints(:,2),'ro')
rectangle('Position',region,'EdgeColor','r','LineWidth',2)
set(gca,'xticklabel',[])
set(gca,'xtick',[])
set(gca,'yticklabel',[])
set(gca,'ytick',[])
xlabel('z/c', 'fontsize',fsl,'FontName','Times-Roman')
ylabel('x/c', 'fontsize',fsl,'FontName','Times-Roman')
title('Dewarped image', 'fontsize',fsl,'FontName','Times-Roman')
% 
% subplot(1,2,2)
% pcolor(Z,X,Area);   
% shading('interp');
% hold on
% cbar = colorbar('eastoutside');
% cbar.Label.String = 'Intensity (counts)';
% cbar.FontSize = fsl;
% %cbar.Ticks = [];
% colormap(gca,jet(128))
% set(gca,'clim',lims)
% axis tight equal
% set(gca,'YDir','reverse')
% set(gca, 'Layer', 'top')
% ylim([x_upstream x_downstream])
% xlim([-0.3/2 0.3/2])
% xlabel('z/c', 'fontsize',fsl,'FontName','Times-Roman')
% ylabel('x/c', 'fontsize',fsl,'FontName','Times-Roman')
% set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
% 
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1)-400 pos(2)-300 width*260, height*120]); %<- Set size

%%
%Spanwise average
Avg = mean(Area,2);
Gradient = gradient(Avg);

figure
yyaxis left
plot(xc,Avg,'LineWidth',2)
ylabel('Intensity (counts)', 'fontsize',fsl,'FontName','Times-Roman')
yyaxis right
plot(xc,Gradient,'LineWidth',2)
grid on
xlabel('x/c', 'fontsize',fsl,'FontName','Times-Roman')
ylabel('Numerical gradient', 'fontsize',fsl,'FontName','Times-Roman')
title('Spanwise-averaged', 'fontsize',fsl,'FontName','Times-Roman')
set(gca, 'FontSize', fsz, 'LineWidth', alw,'Box','off','XMinorTick','on','YMinorTick','on','TickLength',[0.02 0.02]);
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*130, height*120]); %<- Set size

%% saving
save(['Mapped_no_flow_' num2str(Im_number) '.mat']);




