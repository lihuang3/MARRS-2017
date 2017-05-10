clear
clc
load('VeinsInBrain.mat')

figure(4);clf
mapshow(x_obs,y_obs,'DisplayType','polygon',...
'FaceColor',[182,228,255]./255,'LineStyle','none')

polyBin(6,4,6400)=0;
polyBin(1,1,:)=1;
for ii = 1:numel(x_obs)-1
    if (ii ~=74) && (ii+1 ~=74)
        Bx1 = 1+floor(x_obs(ii)/5);
        By1 = 1+floor(y_obs(ii)/5);
        Bx2 = 1+floor(x_obs(ii+1)/5);
        By2 = 1+floor(y_obs(ii+1)/5);
        for jx = min(Bx1,Bx2):max(Bx1,Bx2)
            for jy = min(By1,By2):max(By1,By2)
            
                cnt = polyBin(1,1,jx+80*(jy-1));
                
                polyBin(cnt+1,:,jx+80*(jy-1)) = [x_obs(ii),y_obs(ii),x_obs(ii+1),y_obs(ii+1)];
                
                polyBin(1,1,jx+80*(jy-1)) = cnt +1;
            end
            
        end
        
    end
end