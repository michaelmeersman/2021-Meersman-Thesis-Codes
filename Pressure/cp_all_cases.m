clear; close
width = 4.9;
height = 4.5;

testcase(1) = load('cp_selected_case_1.mat');
testcase(2) = load('cp_selected_case_2.mat');
testcase(3) = load('cp_selected_case_3.mat');
testcase(4) = load('cp_selected_case_4.mat');
testcase(5) = load('cp_selected_case_5.mat');

figure

colors = [0.9290 0.6940 0.1250;
          138, 14, 6;
          0, 24, 102;
          0.4940 0.1840 0.5560;
          48, 171, 0]/255;

for i = [2 3 5]
    plot(testcase(i).upper_mid_taps,testcase(i).cp_upper_mid,'Color',colors(i,:),'LineWidth',1.5,...
        'DisplayName',['Flight, $Re_{c} = $' num2str(round(testcase(i).Re,-3), '%.f') ', $\alpha = $' num2str(testcase(i).alpha_xfoil)])
    hold on
    plot(testcase(i).xfoil_x_c,testcase(i).xfoil_Cp,'-.','Color',colors(i,:),'LineWidth',1.5,...
        'DisplayName',['XFOIL, $Re = $' num2str(testcase(i).Re_r) ', $\alpha = $' num2str(testcase(i).alpha_xfoil)])
    qq = scatter(testcase(i).upper_mid_taps,testcase(i).cp_upper_mid,30,'filled','MarkerFaceColor',colors(i,:));
    qq.Annotation.LegendInformation.IconDisplayStyle = 'off';
    qqq = plot(testcase(i).lower_mid_taps,testcase(i).cp_lower_mid,'Color',colors(i,:),'LineWidth',1.5);
    qqq.Annotation.LegendInformation.IconDisplayStyle = 'off';
    qqqq = scatter(testcase(i).lower_mid_taps,testcase(i).cp_lower_mid,30,'filled','MarkerFaceColor',colors(i,:));
    qqqq.Annotation.LegendInformation.IconDisplayStyle = 'off';
    qqqqq = plot(testcase(i).x_close,testcase(i).cp_close,'LineWidth',1.5,'Color',colors(i,:));
    qqqqq.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     ee = errorbar(testcase(i).upper_mid_taps,testcase(i).cp_upper_mid,testcase(i).err_upper,'k','LineStyle','none');
%     ee.Annotation.LegendInformation.IconDisplayStyle = 'off';
%     eee = errorbar(testcase(i).lower_mid_taps,testcase(i).cp_lower_mid,testcase(i).err_lower,'k','LineStyle','none');
%     eee.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
legend('interpreter','latex','location','south','FontSize',15)
legend show



ylim([min([testcase(5).xfoil_Cp])-0.2 max([testcase(5).xfoil_Cp])])
xlim([-0.05 1])
set(gca, 'YDir','reverse')
set(gca, 'LineWidth', 1.5, 'FontSize', 15)
xlabel('$x/c$','Interpreter','latex','FontSize',24)
ylabel('$C_p$','Interpreter','latex','FontSize',24)
set(gca, 'FontName','Times','FontWeight','bold')
ax.LabelFontSizeMultiplier = 20;
% title('NACA $64_{3}-618$ Flight Test $C_p$','Interpreter','latex','FontSize',24)
grid on
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-500 width*180, height*150]); %<- Set size
hold off