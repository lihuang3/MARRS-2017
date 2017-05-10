
load('ExptData1210.mat');
X = (linspace(1,10,10));
ii = 10;
figure(10);hold on

errorbar(X(1:ii),mean(ExptData(:,1:ii)),std(ExptData(:,1:ii)), 'LineWidth',2)


% 
% load('ExptData1111.mat');
% errorbar(X(1:ii),mean(ExptData(:,1:ii)),std(ExptData(:,1:ii)), 'LineWidth',2)
% 


% xlabel('Number of Particles')
% ylabel('Number of Steps')
% legend('Obstacle-Weighted RRT','Location','northwest')
% %title('Collect the Farthest Particles First')
% 
% % axis([0 1200 0 7000])
% % xticks=[2 128 256 512 1024 ];
% % xticklabels={'2' '128' '256' '512' '1024'};
% 
% set(gca,'FontSize',16,'XTick',xticks,'XTicklabels',xticklabels);

% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 9])
% print -djpeg DataProcess.jpg -r400
%% Y = 2.^(log(X).^0.2)