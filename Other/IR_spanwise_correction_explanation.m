clear; close all

testcase = 5;
load(['corrected_case_' num2str(testcase) '.mat']);



span = Time(30,:);
span2 = double(Mapped.IM_10.Area(30,:));
xc = Mapped.IM_10.xc;
zc = Mapped.IM_10.zc;



third = polyfit(zc,span2,3);
fitted = polyval(third,zc);

fitted_n = fitted-max(fitted,[],2);
span2_corrected = span2 - fitted_n;


h(1) = figure;
subplot(1,7,1:4)
width = 7.5;
height = 4;
plot(zc,span2,'Linewidth',1.5)
hold on
plot(zc,fitted,'Linewidth',1.5)
plot(zc,span2_corrected,'Linewidth',1.5)
xlim([-.155 .155])
set(gca,'Fontsize',14)
ylabel('Intensity Counts','Fontsize',18,'interpreter','latex')
xlabel('z/c','Fontsize',18,'interpreter','latex')
title('Example spanwise correction','Fontsize',18,'interpreter','latex')
legend('Spanwise Distribution','3rd Order Fit','Corrected Distribution',...
    'interpreter','latex','location','southwest','Fontsize',16)



subplot(1,7,5:7)
pcolor(Mapped.IM_10.zc,Mapped.IM_10.xc,Mapped.IM_10.Area)
hold on
set(gca, 'ydir','reverse')
shading('interp');
cbar = colorbar('eastoutside');
cbar.Label.String = 'Intensity';
cbar.FontSize = fsl;
cbar.Label.Interpreter = 'latex';
cbar.Ticks = [];
set(gca,'XTick',[],'YTick',[])
colormap(gca,jet(128))
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size
yline(xc(30),'Linewidth',3)
ax1 = gca;



h(2) = figure;
subplot(1,7,1:4)
plot(X(:,1),No_flow_avg,'Linewidth',1.5,'Displayname','Mean chordwise-distribution')
xlim([0.4 0.84])
set(gca,'Fontsize',14)
legend('interpreter','latex','location','southwest','Fontsize',16)
ylabel('Intensity Counts','Fontsize',18,'interpreter','latex')
xlabel('x/c','Fontsize',18,'interpreter','latex')
title('No-flow chordwise correction offset','Fontsize',18,'interpreter','latex')

subplot(1,7,5:7)
pcolor(Z(1,:),X(:,1),No_flow.Interpolated)
hold on
set(gca, 'ydir','reverse')
shading('interp');
cbar = colorbar('eastoutside');
cbar.Label.String = 'Intensity';
cbar.FontSize = fsl;
cbar.Label.Interpreter = 'latex';
cbar.Ticks = [];
set(gca,'XTick',[],'YTick',[])
colormap(gca,jet(128))
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size
ax2 = gca;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILE SAVING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = ['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\Figures\'...
    'Image Correction Explanation\'];
exportgraphics(h(1),[filepath 'spanwise_correction.png'],'Resolution',1200);
exportgraphics(h(2),[filepath 'no_flow_offset.png'],'Resolution',1200);
exportgraphics(h(1),[filepath 'spanwise_correction.eps']);
exportgraphics(h(2),[filepath 'no_flow_offset.eps']);
savefig(h(1),[filepath 'spanwise_correction.fig']);
savefig(h(2),[filepath 'no_flow_offset.fig']);









