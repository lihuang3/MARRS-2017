
%% Demonstrate some other example plots
%generate generic datasets
% numpts = 50;
% rand_x = 2+randn(numpts ,1);
% rand_y = -0.4+2*randn(numpts ,1);
figure(10);hold on
%  load('ExptData1210.mat')    % Std map
%  load('ExptData_owRRRT1214.mat')  % T-map
 load('ExptVB0103.mat')  % Veins in brain

colorArray = [1.*ones(10,1),0.*ones(10,1),0.*ones(10,1)];
%colorArray = [1 0 0];
xgap = 1;
xoffset = -0.1;
[v1,p1] =violin(ExptData,'facecolor',colorArray,'xgap',xgap,'xoffset',xoffset,'medc','','mc','-k*');
legend off
axis([0 2 0 3000])
% 
 load('ExptDataVB0104.mat')   % Veins in brain
%  load('ExptData1111.mat')    % Std map
%    load('ExptDataowRRT1213.mat')   % T-map
%make a violin plot using the violin FEX function
figure(10); 
%subplot(1,2,1)
xgap = 1;
xoffset = 0;
colorArray = [0.*ones(10,1),1.*ones(10,1),0.*ones(10,1)];
[v2,p2] = violin(ExptData,'facecolor',colorArray,'xgap',xgap,'xoffset',xoffset,'medc','','mc','-k*');
legend off

figure(10)

%  load('ExpData10.mat') % Std map
%  load('ExptData1214TRRT.mat') % T-map
 load('Expt0202.mat') % VB

xgap = 1;
xoffset = 0.1;
%  colorArray = [0.*ones(10,1),0.*ones(10,1),1.*ones(10,1)];
 colorArray = [0.*ones(8,1),0.*ones(8,1),1.*ones(8,1)];

ExptData = ExptData(:,1:8);
[v3,p3] = violin(ExptData,'facecolor',colorArray,'xgap',xgap,'xoffset',xoffset,'medc','','mc','-k*');
legend off

%ylim([0 4000])

% legend([v1(1),v2(1)],{'Recursive Method',...
%     'Heuristic Algorithm'},'Location','NorthWest');

legend([v1(1),v2(1),v3(1)],{'Divide-and-conquer with OWRRT',...
    'Heuristics with OWRRT','Heuristics with RRT'},'Location','NorthWest');

axis([0.5 10.5 -10 15000]);

xlabel('Number of microrobots')
ylabel('Number of Steps')
ax = gca;

ax.XTick=linspace(1,10,10);
ax.XTickLabel={'2^1','2^2','2^3','2^4','2^5','2^6','2^7','2^8','2^9','2^{10}'};