track = [];

[RegionVal, RegionId] = sort(diag(Pathway(BrchPts0(:,1),BrchPts0(:,2))),'Descend');



for ii =1:length(EdPts0)
    temp1 = EdPts0(ii,1:2);
    temp2 = temp1;
    cost1 = Pathway(temp2(1),temp2(2));
    cost2 = cost1;
    cost = cost1;
    cost_diff = 0;
    flag = 0;
    
    while flag<=1 && cost1>100
        temp1 = temp2;
        cost1 = cost2;
        cost_diff = 0;
        cost = cost1;
        if cost1 > 100
            if flag==0
                while(cost_diff<26 && Brch0(temp2(1),temp2(2))==0)
                    nb = Pathway(-1+temp2(1):1+temp2(1),-1+temp2(2):1+temp2(2));
                    nb(2,2)=0;
                    if ~isempty(find(nb==cost)) 
                        [row col] = find(nb==cost);
                        row = row-2;
                        col = col-2;
                        temp = [row(1)+temp2(1),col(1)+temp2(2)];
                        if Brch0(temp(1),temp(2))==0
                           [row col] = find(nb==cost-1);
                           row = row-2;
                           col = col-2;
                           temp2 = [row(1)+temp2(1),col(1)+temp2(2)];
                        else
                           temp2= temp;

                        end
                    else
                        [row col] = find(nb==cost-1);
                        row = row-2;
                        col = col-2;
                        temp2 = [row(1)+temp2(1),col(1)+temp2(2)];
                    end
                    cost2 = Pathway(temp2(1),temp2(2));
                    cost_diff = cost1-cost2;
                    cost = cost-1;
                end
            else
                 while cost_diff<26
                    nb = Pathway(-1+temp2(1):1+temp2(1),-1+temp2(2):1+temp2(2));
                    nb(2,2)=0;
                    if ~isempty(find(nb==cost)) 
                        [row col] = find(nb==cost);
                        row = row-2;
                        col = col-2;
                        temp = [row(1)+temp2(1),col(1)+temp2(2)];
                        if Brch0(temp(1),temp(2))==0
                           [row col] = find(nb==cost-1);
                           row = row-2;
                           col = col-2;
                           temp2 = [row(1)+temp2(1),col(1)+temp2(2)];
                        else
                           temp2= temp;

                        end
                    else
                        [row col] = find(nb==cost-1);
                        row = row-2;
                        col = col-2;
                        temp2 = [row(1)+temp2(1),col(1)+temp2(2)];
                    end
                    cost2 = Pathway(temp2(1),temp2(2));
                    cost_diff = cost1-cost2;
                    cost = cost-1;
                end
                
            end
            
            angle = atan2(double(temp2(1)-temp1(1)),double(temp2(2)-temp1(2)));
            track = [track;angle, double(cost_diff)];
        end
        if Brch0(temp2(1),temp2(2))==1
            flag = flag + 1/2;
        else
            if flag > 0
                flag = flag + 1/2;
            end
        end
    end
    
    
    track = [track;NaN, NaN];

end