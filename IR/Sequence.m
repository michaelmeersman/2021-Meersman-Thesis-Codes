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
setpoint = '20210918_105310_SEQ_A.tiff';
Frames = (100:360); %frames to read
lims = [7100 7900]; %for plotting
ROI = [1 512; 1 640];
%%
%gif
for k=Frames

    Im = imread(setpoint,k);
    Image = Im(ROI(1,1):ROI(1,2),ROI(2,1):ROI(2,2));

    h = figure;
    imagesc(Image)
    colormap(gray)
    cbar = colorbar;
    cbar.Label.String = 'Intensity (counts)';
    cbar.FontSize = fsz;
    cbar.FontName = 'Times-Roman';
    cbar.Location = 'eastoutside';
    set(gca,'clim',lims)
    axis tight equal
    set(gca,'xticklabel',[])
    set(gca,'xtick',[])
    set(gca,'yticklabel',[])
    set(gca,'ytick',[])
    xlabel('z/c', 'fontsize',fsl,'FontName','Times-Roman')
    ylabel('x/c', 'fontsize',fsl,'FontName','Times-Roman')
    title(num2str(k), 'fontsize',fsl,'FontName','Times-Roman')
    set(gcf,'color','w');
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*140, height*140]); %<- Set size
    
    %Phase averaged combined Gifs
    output_raw = 'FT_SEQ_fullimage.gif';
    % Capture the plot as an image 
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if k == Frames(1)
      imwrite(imind,cm,output_raw,'gif', 'Loopcount',inf,'DelayTime',0.4); 
    else 
      imwrite(imind,cm,output_raw,'gif','WriteMode','append','DelayTime',0.4); 
    end 
end