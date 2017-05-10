figure(2);clf
axis([-0.1 100.1 -0.1 100.1]);

x_obs = [0 18 23 0 NaN 28 65 65 33 28 NaN 100 75 75 100 NaN ...
   0 18 18 0 NaN 28 45 33 28 28 NaN 55 65 65 43 55 NaN ...
   0 28 28 0 NaN 100 38 38 100 NaN 0 80 80 0];
y_obs = [90 90 75 75 NaN 90 90 75 75 90 NaN 90 90 40 40 NaN ...
   65 65 40 40 NaN 65 65 40 40 65 NaN 65 65 40 40 65 NaN ...
   30 30 15 15 NaN 30 30 15 15 NaN 5 5 0 0];

    
mapshow(x_obs,y_obs,'DisplayType','polygon',...
    'FaceColor',[182,228,255]./255,'LineStyle','none')

hold on

mapshow([0 100 100 0 0],[0 0 100 100 0],'LineWidth',2,'Color','black')