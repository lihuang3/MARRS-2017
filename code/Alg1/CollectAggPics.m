
% axis([-0.1 400.1 -0.1 314.1])
% axis on
%mapshow([0 220 220 0 0],[0 0 140 140 0],'LineWidth',2,'Color','black')
% xlabel('X')
% ylabel('Y')
%legend('Orig Data','Filtered Data','Location','southeast')
%title('')

% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];

addpath('C:\Users\lhuang28\Documents\GitHub\MagneticController\lihuang\Motion Planning\RRT\Alg1\pic')
fig = gcf;
axis off
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];


set(gca,'FontSize',26)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 9])
Num = 1;
filename = sprintf('Agg%02d.jpg',Num);
while exist(filename)
    Num = Num + 1;
    filename = sprintf('Agg%02d.jpg',Num);
end

print('-djpeg', filename, '-r400');