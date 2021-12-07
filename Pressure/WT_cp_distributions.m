clear; close all

width = 6;
height = 5;

xpos = readmatrix(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    '2.1 All WT Data\xpos.csv']);
WT_x_c = xpos/12;
WT_x_c([49 50]) = [];

top1 = 1:7;
top2 = 8:34;
top3 = 35:41;
bottom1 = 42:46;
bottom2 = 47:53;
bottom3 = 56:60;
static = 61;
total = 62;
k = 0;
for alpha_wt = [0 2 4 6]
    k = k+1;
    WT = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
        '2.1 All WT Data\Taps_AoA_' num2str(alpha_wt) '_Re_200k.CSV']);
    WT_p = WT.data(:,3:62);
    WT_p(:,[49 50]) = [];
    WT_static = WT.data(:,64);
    WT_total = WT.data(:,65);
    WT_q = 0.5*1.11*10.8^2;
    WT_Cp = (mean(WT_p)-mean(WT_static))*6894.76/WT_q;
    
    wing = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
        '2.1 All WT Data\Modified_NACA _643_618.dat']);
    wing = wing.data;
    wing([44 87],:) = [];
    [polar, solution] = xfoil(wing,alpha_wt,200000,0,'ppar n 200','oper/vpar n 10','oper/iter 500');
    xfoil_x_c = solution.xcp;
    xfoil_Cp = solution.cp;
%     x_close = [upper_mid_taps(end) lower_mid_taps(end)];
%     cp_close = [cp_upper_mid(end) cp_lower_mid(end)];
    
    
    p(k) = plot(WT_x_c(horzcat(top2,fliplr(bottom2))),WT_Cp(horzcat(top2,fliplr(bottom2))),...
        '-o','LineWidth',1.5,'DisplayName',['$\alpha = $' num2str(alpha_wt)]);
    hold on
    c = get(p(k),'color');
    pp(k) = plot(xfoil_x_c,xfoil_Cp,'-.','Color',c,'LineWidth',1);
    pp(k).Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    
end

xlim([-0.1 1])
ylim([-1.7 1.2])
set(gca, 'YDir','reverse')
set(gca, 'LineWidth', 1.5, 'FontSize', 15)
legend('location','south','FontSize',20,'interpreter','latex')
xlabel('$x/c$','Interpreter','latex','FontSize',24)
ylabel('$C_p$','Interpreter','latex','FontSize',24)
% set(gca, 'FontName','Times','FontWeight','bold')
ax.LabelFontSizeMultiplier = 20;
title('$Re=200k$','interpreter','latex','FontSize',24)
grid on
hold off
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-500 width*120, height*120]); %<- Set size
