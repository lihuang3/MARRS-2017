[row col] = find(Pathway==max(Pathway(:)));
for jj = 1:numel(row)
   if Ed(row(jj),col(jj))==0
       Ed(row(jj),col(jj))=1;
   end
    
end

Ed0 = Ed;
Brch0 = Brch;
Ed0(goal(1),goal(2)) = 1;
Brch0(goal(1),goal(2)) = 1;

 [row col] = find(Ed==1);
 EdPts = [row col];
 
  [row col] = find(Brch0==1);
 BrchPts0 = [row col];

 BEdist = 30;% 15;
 BBdist = 40; %15;
 
for ii = 1:length(EdPts)
    temp1 = EdPts(ii,1:2); % Starting point
    temp = temp1;
    cost1 = Pathway(temp(1),temp(2));   % Starting point cost
    while(Brch(temp(1),temp(2))==0)
        cost = Pathway(temp(1),temp(2));
        nb = Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
        nb(2,2)=0; 
        if ~isempty(find(nb==cost))  
            [row col] = find(nb==cost);
            row = row-2;
            col = col-2;
            
            if Brch(row(1)+temp(1),col(1)+temp(2))==1
                temp = [row(1)+temp(1),col(1)+temp(2)];
                cost2 = Pathway(temp(1),temp(2));
            else
                [row col] = find(nb==cost-1);
                row = row-2;
                col = col-2;
                temp = [row(1)+temp(1),col(1)+temp(2)];
                cost2 = Pathway(temp(1),temp(2));
            end
        else
            [row col] = find(nb==cost-1);
            row = row-2;
            col = col-2;
            temp = [row(1)+temp(1),col(1)+temp(2)];
            cost2 = Pathway(temp(1),temp(2));
        end
    end
   if (cost1-cost2<BEdist)

        Brch0(temp(1),temp(2))=0;
        Ed0(temp(1),temp(2))=1;
        Ed0(temp1(1),temp1(2))=0;
   end
       

end

[row col] = find(Ed0);
EdPts0 = [row col];

for ii = 1:size(EdPts0,1)
    tempEd = EdPts0(ii,:);
    if Ed0(tempEd(1),tempEd(2))==0
        continue
    end
    temp = tempEd;
    cost = Pathway(temp(1),temp(2));
    while 1     
        nb = Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
        nb(2,2) = 0;
        if ~isempty( find(nb==cost))
            [row col] = find(nb==cost);

            row = row-2;
            col = col-2;
            
            
            
            if Ed0(temp(1)+row(1),temp(2)+col(1))==1
                
                temp = [temp(1)+row(1),temp(2)+col(1)];
                Ed0(temp(1),temp(2)) = 0;
            else
                [row col] = find(nb==cost-1);
                if isempty(row)
                   break 
                else
                    row = row-2;
                    col = col-2;
                    temp = [row(1)+temp(1),col(1)+temp(2)];
                    cost = cost - 1;
                    if Ed0(temp(1),temp(2))==1
                        Ed0(temp(1),temp(2)) = 0;
                    end
                end
            end
        else
            [row col] = find(nb==cost-1);
            if isempty(row)
                break
            else
                row = row-2;
                col = col-2;

                temp = [row(1)+temp(1),col(1)+temp(2)];
                cost = cost - 1;
                if  Ed0(temp(1),temp(2))==1
                    Ed0(temp(1),temp(2)) = 0;
                    break
                end
            end
        end
    end 
end

%% Combine Brch
for ii = 1:size(BrchPts0,1)
    tempBh = BrchPts0(ii,:);
    if Brch0(tempBh(1),tempBh(2))==0
        continue
    end
    temp = tempBh;
    cost1 = Pathway(temp(1),temp(2));
    cost2 = cost1;
    cost = cost2-cost1;
    while cost<BBdist && cost2>100    
        nb = Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
        nb(2,2) = 0;
        if ~isempty( find(nb==cost2))
            [row col] = find(nb ==cost2);
            row = row-2;
            col = col-2;
            if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                temp = [temp(1)+row(1),temp(2)+col(1)];
                Brch0(tempBh(1),tempBh(2))=0;
                break
            else
                [row col] = find(nb ==cost2-1);
                row = row-2;
                col = col-2;
                cost2 = cost2-1;
                cost = cost1-cost2;
                if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                    temp = [temp(1)+row(1),temp(2)+col(1)];
                    Brch0(tempBh(1),tempBh(2))=0;
                    break
                else   

                   temp =[temp(1)+row(1),temp(2)+col(1)]; 
                end
                
            end
        else
            [row col] = find(nb ==cost2-1);
            row = row-2;
            col = col-2;
            cost2 = cost2-1;
            cost = cost1-cost2;
            if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                temp = [temp(1)+row(1),temp(2)+col(1)];
                Brch0(tempBh(1),tempBh(2))=0;
                break
            else   
               
               temp =[temp(1)+row(1),temp(2)+col(1)]; 
            end
            
            
        end
        
        
    end
end


Brch0(goal(1),goal(2))=1;
[row col] = find(Brch0);
BrchPts0 = [row col];
Ed0(goal(1),goal(2)) = 1;
[row col] = find(Ed0);
EdPts0 = [row col];


figure
imshow(Skel);
hold on
%plot(BrchPts(:,2),BrchPts(:,1),'r.');
f1 = scatter(BrchPts0(:,2),BrchPts0(:,1),20,'filled');
f1.CData = [1 0 0];
f2 = scatter(EdPts0(:,2),EdPts0(:,1),20,'filled');
f2.CData = [0 0 1];
if ~isempty(sdEdPts)
    f3 = scatter(sdEdPts(:,2),sdEdPts(:,1),5,'filled');
    f3.CData = [0,0,0.8];
end
