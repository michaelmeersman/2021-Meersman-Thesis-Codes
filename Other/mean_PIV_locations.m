clear; close

width = 5;
height = 4;

data = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
                  'Extracted Data\PIV Mean Bubble Locations.csv'],',',2);
              
sep_x = data.data(:,3);              
sep = data.data(:,4);  

trans_x = data.data(:,5);
trans = data.data(:,6);

reatt_x = data.data(:,1);
reatt = data.data(:,2);


p1 = plot(sep_x,sep,':','Color',[0, 119, 237]/256,'LineWidth',1.5);
hold on
p2 = plot(trans_x,trans,':r','LineWidth',1.5);
p3 = plot(reatt_x,reatt,':','Color',[252, 186, 3]/256,'LineWidth',1.5);
scatter(sep_x,sep,'filled','MarkerFaceColor',[0, 119, 237]/256,'LineWidth',1.5,'MarkerEdgeColor','k')
scatter(trans_x,trans,'r','filled','LineWidth',1.5,'MarkerEdgeColor','k')
scatter(reatt_x,reatt,'filled','MarkerFaceColor',[252, 186, 3]/256,'LineWidth',1.5,'MarkerEdgeColor','k')
set(gca,'FontSize',14)
xlabel('$\alpha\:^{\circ}$','interpreter','latex','FontSize',20)
ylabel('$x/c$','interpreter','latex','FontSize',20)

p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
p3.Annotation.LegendInformation.IconDisplayStyle = 'off';

legend('Separation','Transition','Reattachment','FontSize',16,'interpreter','latex')

grid on

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size
