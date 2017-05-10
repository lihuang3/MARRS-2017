ArraySize = sum(Skel(:));
Pathway = uint16(Skel);
cost = 1e2;
%Pathway(245,273)=cost;
goal = [331,333]; % StdMap
%goal=[297,169];  % VB
%goal = [177,244];

%goal=[280,327];  % T
Pathway(goal(1),goal(2)) =cost;
frtr = goal;
%frtr = [245,273];
rg = 1;
ii = 1;
jj = 1;
sdEdPts = [];
while ~isempty(frtr) %(jj<ArraySize)
    temp = Pathway(-1+frtr(ii,1):1+frtr(ii,1),-1+frtr(ii,2):1+frtr(ii,2));
    [row col] = find(temp==1);
    row = row-2;
    col = col-2;
    cnt = numel(row);
    if cnt>0
        Pathway(-1+frtr(ii,1):1+frtr(ii,1),-1+frtr(ii,2):1+frtr(ii,2))=...
            -Pathway(-1+frtr(ii,1):1+frtr(ii,1),-1+frtr(ii,2):1+frtr(ii,2)).*uint16(blkdiag(0,1,0))+...
            uint16(temp ==1).*Pathway(frtr(ii,1),frtr(ii,2))...
            +Pathway(-1+frtr(ii,1):1+frtr(ii,1),-1+frtr(ii,2):1+frtr(ii,2));

        cost = cost +1;
        for ij = 1:cnt
            frtr = [frtr; row(ij)+frtr(ii,1), col(ij)+frtr(ii,2)];     
        end
    else
       if Ed(frtr(ii,1),frtr(ii,2))==0
            sdEdPts = [sdEdPts;frtr(ii,1:2)]; 
       end
    end
    frtr(ii,:) = [];
    jj = jj + cnt;
end
if ~isempty(sdEdPts)
    sdEdPts = sdEdPts(find(diag(Brch(sdEdPts(:,1),sdEdPts(:,2)))==0),:);
end