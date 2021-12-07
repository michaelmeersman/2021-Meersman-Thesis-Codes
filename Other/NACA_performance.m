clear; close


NACA = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
    'Extracted Data\NACA_6_series_performance.csv'],',');

alpha = NACA.data(:,1);
Cl = NACA.data(:,2);
Cd = NACA.data(:,3);


width = 4;
height =3;


plot(Cl,Cd,'k','Linewidth',1.5)
%\alpha\:^{\circ}
xlabel('$C_l$','Fontsize',20,'interpreter','latex')
ylabel('$C_d$','Fontsize',20,'interpreter','latex')
title('NACA $64_{3}-618$, Re=500k, Ncrit = 5','Fontsize',20,'interpreter','latex')
grid on
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size