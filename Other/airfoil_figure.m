clear; close

AF = readmatrix(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\'...
    '4. Airfoils for Solidworks\Modified 6 Series Airfoil Coordinates (checked 2_23_21).xlsx'],...
    'range','A1:B86');

x = AF(:,1);
y = AF(:,2);

width = 6;
height =3;


plot(x,y,'k','Linewidth',1.5)
hold on
plot([x(end) x(1)],[y(end) y(1)],'k','Linewidth',1.5)
axis equal
ylim([-0.2 0.2])
xlim([0 1])
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-300 width*120, height*120]); %<- Set size