clear; close

% Endevcos are hardware low-pass filtered at 10k
for AOA = 0:10
    for RE = 1:2
        Re = (RE+1)*100;
        E(RE,AOA+1) = load(['Endevco_RE_' num2str(Re) 'k_AOA_' num2str(AOA) '.mat']);
    end
end

%%
close all
pos = [0.35 0.5 0.65];
u_inf = [10.7 16.2]; %m/s
c = 0.305; %m
str = [16 20];

x = linspace(1,10000);
y = ((x).^-5/3)*15000;

for q = 1:2
    figure(q)
    for i = [0 2 4 6]+1
        for j = 3 % selection of endevco sensor
            loglog(E(q,i).power_freqs(j,:)*c/u_inf(q),E(q,i).power_mags(j,:),'DisplayName',...
                ['$\alpha = $ ' num2str(i-1)],...
                'LineWidth',1.4)
            hold on
        end
    end
    % Mark the "mean" strouhal number
    xx = xline(str(q),'k','Label',['$St\approx$' num2str(str(q))],...%'Color',[199, 0, 0]/256,...
        'LineWidth',2,'interpreter','latex','LabelVerticalAlignment','top','FontSize',18);
    xx.Annotation.LegendInformation.IconDisplayStyle = 'off';
    xxx = loglog(x,y,'Linewidth',1.2,'Color','k','DisplayName','$f^{-5/3}$');
%     xxx.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    width = 5.5;
    height = 4.5;
    
    % legend('0.35c','0.5c','0.65c','Interpreter','latex')
    legend('Location', 'northwest','FontSize',14,'Interpreter','latex','Numcolumns',2)
    xlab = xlabel('$\frac{f \cdot c}{u_\infty}$','Interpreter','latex','FontSize',24);
    ylab = ylabel('PSD','Interpreter','latex','FontSize',24);
    xlab_pos = get(xlab,'Position');
%     set(xlab,'Position',[xlab_pos(1) xlab_pos(2) xlab_pos(3)]);
    % xx = xline(500,'-.k','LineWidth',1.2);
    % xx.Annotation.LegendInformation.IconDisplayStyle = 'off';
    % xxx = xline(625,'-.k','LineWidth',1.2);
    % xxx.Annotation.LegendInformation.IconDisplayStyle = 'off';
    left = 700;
    right = 1500;
    x_loc = [left right right left];
    y_loc = [10 10 10^-25 10^-25];
    % ppp = patch(x_loc,y_loc,'r','facealpha',0.1);
    % ppp.Annotation.LegendInformation.IconDisplayStyle = 'off';
    xlim([3 50])
    ylim([4e-9 2e-3])
    set(gca,'FontSize',18)
    set(xlab,'FontSize',28)
    set(ylab,'FontSize',20)
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size
    grid on
    
    
end




