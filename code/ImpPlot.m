track0 = track;
[RegionVal, RegionID] = sort(diag(Pathway(BrchPts0(:,1),BrchPts0(:,2))),'Descend');

ii = 1;
jj = 6;
%jj = 6;
while ~isnan( track0(ii,2.*RegionID(jj)-1))
    ii = ii +1; 
end

tar = rtrack(1,2.*RegionID(jj)-1);
tar = [tar,  sum(track0(ii-2:ii-1,2*RegionID(jj)))];
tar = [tar, rtrack(2,2.*RegionID(jj))];
jj =jj+1;




kk = 1;
cur_track = zeros(1,2*ConnMat(RegionID(jj),1));
ii = 1;
ll = 1;
while ii<= size(track0,1) && (track0(ii,2.*RegionID(jj))>0 || isnan( track0(ii,2.*RegionID(jj)-1)) )
    if ~isnan( track0(ii,2.*RegionID(jj)-1))
        cur_track(ll,2*kk-1:2*kk) = [track0(ii,2.*RegionID(jj)-1),track0(ii,2.*RegionID(jj))];
        ll = ll + 1;
    else
        kk = kk + 1;
        ll = 1;
    end
    ii = ii + 1;
end
    
YC = tar(2).*ones(factorial(ConnMat(RegionID(jj),1)),1);
Y = 0.*YC;


seq = perms(linspace(1,ConnMat(RegionID(jj),1),ConnMat(RegionID(jj),1)));
for ii = 1:factorial(ConnMat(RegionID(jj),1))
    mm = 2;
    for kk = 1:ConnMat(RegionID(jj),1)
        
       ll = 1;
       kk0 = seq(ii,kk);
       while ll<=size(cur_track,1) && cur_track(ll,2*kk0) > 0
           temp = YC(ii,mm-1)+cur_track(ll,2*kk0).*cos(cur_track(ll,2*kk0-1)-tar(1));
           if temp >tar(3)
                temp = tar(3);
           end
           Y(ii,mm) = temp-YC(ii,mm-1);
           YC(ii,mm) = temp;
           mm = mm + 1;
           ll = ll + 1;
       end
    end
end

[~, opt] = max(min(YC,[],2));

% while ii<= size(track0,1) && (track0(ii,2.*RegionID(jj))>0 || isnan( track0(ii,2.*RegionID(jj)-1)) )
%     if ~isnan( track0(ii,2.*RegionID(jj)-1))
%         temp = YC(end)+track0(ii,2.*RegionID(jj)).*cos(track0(ii,2.*RegionID(jj)-1)-tar(1));
%         if temp >tar(3)
%             temp = tar(3);
%         end
%         YC = [YC, temp];
%         Y = [Y,track0(ii,2.*RegionID(jj)).*cos(track0(ii,2.*RegionID(jj)-1)-tar(1))];
%     end
%     ii = ii + 1;
% end
figure
X = linspace(1,size(Y,2),size(Y,2));
stairs(X,Y(opt,:),'LineWidth',2); hold on
plot(X,X.*0, '-r');
xlabel('step number')
ylabel('step length * cos(angle difference)')

figure
XC = X;
stairs(XC,YC(opt,:),'LineWidth',2);hold on
plot(X,X.*0, '-r');
xlabel('step number')
ylabel('accumulated effect')