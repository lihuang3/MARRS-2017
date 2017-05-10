
track0 = track;

[RegionVal, RegionID] = sort(diag(Pathway(BrchPts0(:,1),BrchPts0(:,2))),'Descend');


for ll = 1:size(track0,2)/2
fig1 = figure;
fig1.Position = [200 200 400 400];
imshow(EBW); hold on
trackX = ConnMat(RegionID(ll),3);
trackY = ConnMat(RegionID(ll),2);
ii=1;
while ii<=size(track0,1) && ~(track0(ii,2*RegionID(ll))==0)
    if ~isnan(track0(ii,2*RegionID(ll)))
        trackX = [trackX,trackX(end)+track0(ii,2*RegionID(ll))*cos(track0(ii,2*RegionID(ll)-1))];
        trackY = [trackY,trackY(end)+track0(ii,2*RegionID(ll))*sin(track0(ii,2*RegionID(ll)-1))] ;
    else
        scatter(trackX(end),trackY(end),35,'filled');
        
    end
    ii = ii +1;
end



plot(trackX,trackY,'k','LineWidth',2);
scatter(trackX(1),trackY(1),50)
hold off
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'top';

fig2 = figure;
fig2.Position = [600 200 400 400];
imshow(EBW); 
hold on
trackX = ConnMat(RegionID(ll),3);
trackY = ConnMat(RegionID(ll),2);
w = 10;
w0 = 0;
wL = 17.5;
wR  = 17.5;
markpts = [];
ii = 1;
while ii<=size(track0,1) && ~(track0(ii,2*RegionID(ll))==0)
    if ~isnan(track0(ii,RegionID(ll)*2))
        if isempty(markpts) 
            trackX = [trackX,trackX(end)+track0(ii,2*RegionID(ll))*cos(track0(ii,2*RegionID(ll)-1))];
            trackY = [trackY,trackY(end)+track0(ii,2*RegionID(ll))*sin(track0(ii,2*RegionID(ll)-1))];
        else
            u = track0(ii,2*RegionID(ll)).*[cos(track0(ii,2*RegionID(ll)-1)),sin(track0(ii,2*RegionID(ll)-1))];

           if norm(OrtVecL*u'+w0) <= w/2
               X = trackX(end)+u(1);
               Y = trackY(end)+u(2);
               w0 = w0 + (OrtVecL*u');
           else
                temp =  u.*(max(min(w0 + (OrtVecL*u'),w/2),-w/2)-w0)./(w0+OrtVecL*u');
                w0 = max(min(w0 + (OrtVecL*u'),w/2),-w/2);
                X = trackX(end)+temp(1);
                Y = trackY(end)+temp(2);

           end


           trackX = [trackX,X];
           trackY = [trackY,Y];
            
        end
    else
        scatter(trackX(end),trackY(end),35,'filled');
        if isempty(markpts)
            markpts = mean(track0(ii-2:ii-1,2*RegionID(ll)-1:2*RegionID(ll)));
            OrtVecL = [-sin(markpts(1)),cos(markpts(1))];
            OrtVecR = -OrtVecL;
        end
        
    end
ii=ii+1;
end

plot(trackX,trackY,'k','LineWidth',2);
scatter(trackX(1),trackY(1),50);
hold off
%axis([-10 300 -10 150])
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'top';
disp(' ')
end
