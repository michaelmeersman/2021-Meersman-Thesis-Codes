clear; close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
width = 9;
height = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

taps = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\' ...
    '7. Matlab Pressure Tap Codes\6series_pressure_tap_FLIGHT.csv']);
upper_mid_taps = taps(1:30,1);
lower_mid_taps = taps(1:8,2);

%%% Define airfoil top and bottom indecies
top1 = 1:7;
top2 = 8:34;
top2_plot = 8:34;
top3 = 35:41;
bottom1 = 42:46;
bottom2 = 47:55;
bottom3 = 56:60;

%%% pressure channels
static = 61;
total = 62;


for AOA = -5:15
    if AOA < 0
        aoastr = ['n' num2str(abs(AOA))];
    elseif AOA >= 0
        aoastr = num2str(AOA);
    end
    
    data_300k = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA'...
        '\Summer WT Data 300k\07_06_300k_' aoastr '.txt'],'\t',4);
    WT_Cp_300 = data_300k.data(:,8);
    WT_x_c_300 = data_300k.data(:,3)/12;
    q = data_300k.data(total,8) - data_300k.data(static,8);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%% XFLR %%%%%%%%%%%%%%%%%%%%%%%%%%
    wing = importdata(['C:\Users\micha\Documents\1.RESEARCH\1.FLIGHT\6. DATA\'...
        '2.1 All WT Data\Modified_NACA _643_618.dat']);
    wing = wing.data;
    wing([44 87],:) = [];
    [polar, solution] = xfoil(wing,AOA,300000,0,'ppar n 200','oper/vpar n 7','oper/iter 500');
    xfoil_x_c = solution.xcp;
    xfoil_Cp = solution.cp;
    % x_close = [upper_mid_taps(end) lower_mid_taps(end)];
    % cp_close = [cp_upper_mid(end) cp_lower_mid(end)];
    
    figure
    plot(xfoil_x_c,xfoil_Cp,'-.b','LineWidth',1.5)
    hold on
    plot(WT_x_c_300(horzcat(top2,fliplr(bottom2))),WT_Cp_300(horzcat(top2,fliplr(bottom2))),'--r','LineWidth',1.5)
    legend('XFOIL','Wind Tunnel')
    xlim([-0.05 1])
    set(gca, 'YDir','reverse')
    set(gca, 'LineWidth', 1.5, 'FontSize', 15)
    title(['Cp Comparisson AOA = ' aoastr])
    xlabel('Chordwise position, $x/c$','Interpreter','latex','FontSize',20)
    ylabel('Pressure Coefficient, $C_p$','Interpreter','latex','FontSize',20)
    set(gca, 'FontName','Times','FontWeight','bold')
    ax.LabelFontSizeMultiplier = 20;
    grid on
    pos = get(gcf, 'Position');
    set(gcf, 'Position', [pos(1)-200 pos(2)-500 width*120, height*150]); %<- Set size
    hold off
end


