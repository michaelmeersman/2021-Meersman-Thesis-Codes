clear; close

width = 6;
height = 5;
fsz = 18;


data = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\Extracted Data\'...
                   'Spalart & Strelets Cf and St correlation DATA.csv']);
               
data2 = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\Extracted Data\'...
                   'displacement_momentum_thickness.csv']);


displ = rmmissing(data2.data(:,2));
displ_x = rmmissing(data2.data(:,1));
mom = rmmissing(data2.data(:,4));
mom_x = rmmissing(data2.data(:,3));

x_sm = 0:0.01:displ_x(end);
displ_I = interp1(displ_x,displ,x_sm);
mom_I = interp1(mom_x,mom,x_sm);

% Calculate shape factor at each location
H = displ_I./mom_I;

[tt, trans] = max(H);
TRANS1 = x_sm(trans);

                         
St = rmmissing(data.data(:,4));
St_x = rmmissing(data.data(:,3));
Cf = rmmissing(data.data(:,2));
Cf_x = rmmissing(data.data(:,1));

x_sm2 = St_x(1):0.01:St_x(end);
ii = [St_x St];
ii = ii(28:80,:);
St_I = interp1(ii(:,1),ii(:,2),x_sm2);
St_grad = gradient(St_I);
[stmin, sti] = min(St_grad);
SEP2 = x_sm2(sti);
[stmax, stmi] = max(St_grad);
TRANS2 = x_sm2(stmi);
[stmaxv, maxvi] = max(St_I);
REATT2 = x_sm2(maxvi);


% Find zeros of Cf
k=1;
for i = 1:length(Cf)
    if abs(Cf(i)) < 0.00013
        zeros(k) = i;
        k = k+1;
    end
end

SEP1 = Cf_x(zeros(1));
REATT1 = Cf_x(zeros(end));

%colors
c2 = [196 0 0]/256;
c1 = [0 182 232]/256;
c3 = [153, 140, 3]/256;


figure
yyaxis left
plot(St_x,St,'Color',[0 63 163]/256,'LineWidth',1.5)
hold on
plot(Cf_x,Cf,'-k','LineWidth',1.5)
set(gca, 'FontName','Times','FontWeight','bold','FontSize',fsz-6)
xlabel('x','interpreter','latex','FontSize',fsz+2)
xx = xline(TRANS1,'--','color',c1,'LineWidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','interpreter','latex');
xx.Annotation.LegendInformation.IconDisplayStyle = 'off';
xx2 = xline(TRANS2,'--','color',c1,'LineWidth',1.5,'LabelHorizontalAlignment','left','LabelVerticalAlignment','middle','interpreter','latex');
xx2.Annotation.LegendInformation.IconDisplayStyle = 'off';
xxx = xline(SEP1,'--','color',c2,'LineWidth',1.5,'LabelVerticalAlignment','middle','interpreter','latex');
xxx.Annotation.LegendInformation.IconDisplayStyle = 'off';
xxx2 = xline(SEP2,'--','color',c2,'LineWidth',1.5,'LabelVerticalAlignment','middle','interpreter','latex');
xxx2.Annotation.LegendInformation.IconDisplayStyle = 'off';
xxxx = xline(REATT1,'--','color',c3,'LineWidth',1.5,'LabelVerticalAlignment','middle','interpreter','latex');
xxxx.Annotation.LegendInformation.IconDisplayStyle = 'off';
xxxx2 = xline(REATT2,'--','color',c3,'LineWidth',1.5,'LabelVerticalAlignment','middle','interpreter','latex');
xxxx2.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylll = yline(0,'--k');
ylll.Annotation.LegendInformation.IconDisplayStyle = 'off';
sc1 = scatter(SEP2,St_I(sti),75,'k','LineWidth',1.2);
sc2 = scatter(SEP1,Cf(zeros(1)),75,'k','LineWidth',1.2);
sc3 = scatter(TRANS2,St_I(stmi),75,'^k','LineWidth',1.2);
sc4 = scatter(REATT2,St_I(maxvi),100,'sk','LineWidth',1.2);
sc5 = scatter(REATT1,Cf(zeros(end)),100,'sk','LineWidth',1.2);
sc1.Annotation.LegendInformation.IconDisplayStyle = 'off';
sc2.Annotation.LegendInformation.IconDisplayStyle = 'off';
sc3.Annotation.LegendInformation.IconDisplayStyle = 'off';
sc4.Annotation.LegendInformation.IconDisplayStyle = 'off';
sc5.Annotation.LegendInformation.IconDisplayStyle = 'off';
xlim([0 St_x(end)])
yyaxis right
plot(x_sm,H,'LineWidth',1.5)
sc6 = scatter(TRANS1,tt,75,'^k','LineWidth',1.2);
sc6.Annotation.LegendInformation.IconDisplayStyle = 'off';
ylabel('Shape Factor, H','interpreter','latex','FontSize',fsz)
legend('$2St_{DNS}$','$C_{f DNS}$','$H_{DNS}$','interpreter','latex','FontSize',fsz-4)
% HH.Annotation.LegendInformation.IconDisplayStyle = 'off';

grid on
box on

ax = gca;
ax.YRuler.Exponent = 0;

pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1)-200 pos(2)-500 width*120, height*120]); %<- Set size