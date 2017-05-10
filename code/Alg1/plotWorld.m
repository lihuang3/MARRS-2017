%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotWorld
%%   plot obstacles and path
function plotWorld(world,tree)
  % the first element is the north coordinate
  % the second element is the south coordinate

  N = 10;
  th = 0:2*pi/N:2*pi;
  %figure(1), clf
  figure(1), hold on

  axis([world.SWcorner(1),world.NEcorner(1),...
      world.SWcorner(2), world.NEcorner(2)]);
  %fig1.Position = [-1300 450 400 300];
  %hold on

%% Tree Node & Branches   
%       for ii = 2:size(tree,1)
%         
%         XT = [tree(ii,1),tree((tree(ii,5)),1)];
%         YT = [tree(ii,2),tree((tree(ii,5)),2)];
%         plot(XT,YT,'-k');hold on
%       end
%%  MainStreet
      MajorParent = unique(tree(:,5));
      MainSt = [MajorParent,histc(tree(:,5),MajorParent)];
      MajorParent = MainSt(MainSt(:,2)>=2,1);
      %scatter( tree(MajorParent,1),tree(MajorParent,2),30,'filled');
      MajorParentFilt = MajorParent(tree(MajorParent(:,1),3)<=5);
      
      a = find(abs(tree(:,1)-20.46)<0.01 & abs(tree(:,2)-99.83)<0.01)
      a = find(MajorParentFilt==a)
      MajorParentFilt(a)=[];
      
      h = scatter( tree(MajorParentFilt,1),tree(MajorParentFilt,2),50,'filled'); 
      h.CData = [255 165 0]./255;  
%% Smooth Path
%       MajorParent = unique(tree(:,5));
%       MainSt = [MajorParent,histc(tree(:,5),MajorParent)];
%       MajorParent = MainSt(MainSt(:,2)>=100,1);
%       scatter( tree(MajorParent,1),tree(MajorParent,2),30,'filled');
% 


   mapshow(world.ce(:),world.cn(:),'DisplayType','polygon','FaceColor',[182,228,255]./255,'LineStyle','none')
  

hold on

mapshow([0 100 100 0 0],[0 0 100 100 0],'LineWidth',2,'Color','black')
  targetH = scatter(99,5,50,'filled');
  targetH.CData = [1 0 0];
  % This line was for circlular obstacles
  %mapshow(x_obs,y_obs,'DisplayType','point','Marker','o')
%   X = path(:,2);
%   Y = path(:,1);

%   plot(X,Y,'r','LineWidth',2);
%   
%  set(gca,'linewidth',2)
%  set(gcf,'PaperUnits','inches','PaperPosition',[0 0 9 9])
%  print -djpeg RRT.jpg -r200