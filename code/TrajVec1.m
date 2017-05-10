track = zeros(1,2*size(BrchPts0,1));
trackL = track;
rtrack = track;
MaxStep = 10;

[RegionVal, RegionID] = sort(diag(Pathway(BrchPts0(:,1),BrchPts0(:,2))),'Descend');

ConnMat = zeros(size(BrchPts0,1),11);

for ii = 1: numel(RegionID)
    % Initialization
    cost1 = RegionVal(ii);
    temp1 = BrchPts0(ii,:);
    temp2 = temp1;
    temp = temp1;
    cost2 = cost1;
    cost = 0;
    
    % search for end point of the current junction
    while 1
       nb = Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
       nb(2,2) = 0;
       cost2 = Pathway(temp(1),temp(2));
       if ~isempty(find(nb == cost2))
           [row col] = find(nb == cost2);
           row = row-2;
           col = col-2;
           if Brch0(temp(1)+row(1),temp(2)+col(1))==1
               temp = [temp(1)+row(1),temp(2)+col(1)];
                ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                ConnMat(ID,1) = ConnMat(ID,1)+1;
                ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                break;
           else
              [row col] = find(nb == cost2-1);
              row = row-2;
              col = col-2;
              if isempty(row)
                  break;
              else
                if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                    temp = [temp(1)+row(1),temp(2)+col(1)];
                    ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                    ConnMat(ID,1) = ConnMat(ID,1)+1;
                    ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                    break;
                else
                    temp = [temp(1)+row(1),temp(2)+col(1)];
                end
              end
           end
       else
          [row col] = find(nb == cost2-1);
          row = row-2;
          col = col-2;
          if isempty(row)
              break;
          else
            if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                temp = [temp(1)+row(1),temp(2)+col(1)];
                ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                ConnMat(ID,1) = ConnMat(ID,1)+1;
                ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                break;
            else
                temp = [temp(1)+row(1),temp(2)+col(1)]; 
            end
          end
       end
    end
    
    
    
end

[RegionVal, RegionID] = sort(diag(Pathway(EdPts0(:,1),EdPts0(:,2))),'Descend');

for ii = 1: numel(RegionID)
    % Initialization
    cost1 = RegionVal(ii);
    temp1 = EdPts0(ii,:);
    temp2 = temp1;
    temp = temp1;
    cost2 = cost1;
    cost = 0;
    
    % search for end point of the current junction
    while 1
       nb = Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
              cost2 = Pathway(temp(1),temp(2));

       nb(2,2) = 0;
       if ~isempty(find(nb == cost2))
           [row col] = find(nb==cost2);
           row = row-2;
           col = col-2;
           if Brch0(temp(1)+row(1),temp(2)+col(1))==1 
                temp = [temp(1)+row(1),temp(2)+col(1)];
                ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                ConnMat(ID,1) = ConnMat(ID,1)+1;
                ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                break;
           else
              [row col] = find(nb == cost2-1);
              row = row-2;
              col = col-2;
              if isempty(row)
                  break;
              else
                if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                    temp = [temp(1)+row(1),temp(2)+col(1)];
                    ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                    ConnMat(ID,1) = ConnMat(ID,1)+1;
                    ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                    break;
                else
                    temp = [temp(1)+row(1),temp(2)+col(1)];
                end
              end
           end
       else
          [row col] = find(nb == cost2-1);
          row = row-2;
          col = col-2;
          if isempty(row)
              break;
          else
            if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                temp = [temp(1)+row(1),temp(2)+col(1)];
                ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                ConnMat(ID,1) = ConnMat(ID,1)+1;
                ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                break;
            else
                temp = [temp(1)+row(1),temp(2)+col(1)]; 
            end
          end
       end
    end
    
    
    
end

if ~isempty(sdEdPts)

[RegionVal, RegionID] = sort(diag(Pathway(sdEdPts(:,1),sdEdPts(:,2))),'Descend');

for ii = 1: numel(RegionID)
    % Initialization
    cost1 = RegionVal(ii);
    temp1 = sdEdPts(ii,:);
    temp2 = temp1;
    temp = temp1;
    cost2 = cost1;
    cost = 0;
    
    % search for end point of the current junction
    while 1
       nb = Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
              cost2 = Pathway(temp(1),temp(2));

       nb(2,2) = 0;
       if ~isempty(find(nb == cost2))
           [row col] = find(nb==cost2);
           row = row-2;
           col = col-2;
           if Brch0(temp(1)+row(1),temp(2)+col(1))==1 
                temp = [temp(1)+row(1),temp(2)+col(1)];
                ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                ConnMat(ID,1) = ConnMat(ID,1)+1;
                ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                break;
           else
              [row col] = find(nb == cost2-1);
              row = row-2;
              col = col-2;
              if isempty(row)
                  break;
              else
                if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                    temp = [temp(1)+row(1),temp(2)+col(1)];
                    ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                    ConnMat(ID,1) = ConnMat(ID,1)+1;
                    ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                    break;
                else
                    temp = [temp(1)+row(1),temp(2)+col(1)];
                end
              end
           end
       else
          [row col] = find(nb == cost2-1);
          row = row-2;
          col = col-2;
          if isempty(row)
              break;
          else
            if Brch0(temp(1)+row(1),temp(2)+col(1))==1
                temp = [temp(1)+row(1),temp(2)+col(1)];
                ID = find(BrchPts0(:,1)==temp(1) & BrchPts0(:,2)==temp(2));
                ConnMat(ID,1) = ConnMat(ID,1)+1;
                ConnMat(ID,2*ConnMat(ID,1):2*ConnMat(ID,1)+1) = temp1;
                break;
            else
                
                temp = [temp(1)+row(1),temp(2)+col(1)]; 
            end
          end
       end
    end
    
    
    
end

end
%%
%%
%%

Init = 1;

for jj = 1:size(ConnMat,1)
    kk = 1;
    for ii =1:ConnMat(jj,1)
        temp1 = ConnMat(jj,2*ii:2*ii+1);
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
                    while(cost_diff == 0 || cost_diff<MaxStep && Brch0(temp2(1),temp2(2))==0)
                        
                        nb = Pathway(-1+temp2(1):1+temp2(1),-1+temp2(2):1+temp2(2));
                        
                        nb(2,2)=0;
                        cost = Pathway(temp2(1),temp2(2));
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
                     while cost_diff<MaxStep
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
                %angle = atan2(double(temp2(2)-temp1(2)),double(temp2(1)-temp1(1)));

                if cost_diff > 0
                    
                    track(kk,2*jj-1:2*jj) = [angle, double(cost_diff)];
                    kk = kk+1;
                end
            end
            if Brch0(temp2(1),temp2(2))==1
                flag = flag + 1/2;
            else
                if flag > 0
                    flag = flag + 1/2;
                end
            end
        end


        track(kk,2*jj-1:2*jj) = [NaN, NaN];
        kk = kk+1;

    end
    
end

for ii = 1:size(trackL,2)./2
   
   temp1 = BrchPts0(ii,:);
   if (temp1 == goal) ~= [1 1]  
       row = [];
       jj = 0;
       while isempty(row)
           jj = jj +1;
           row = find(ConnMat(:,2*jj)==temp1(1) & ConnMat(:,2*jj+1)==temp1(2));
       end

       cnt = 1;
       kk = 1;
       while cnt ~= jj
           if  isnan(track(kk,2*row))
               cnt = cnt +1;
           end
           kk = kk +1;
       end
       k1 = kk;
       while ~isnan(track(kk,2*row))
           kk = kk +1;
       end

       k2 = kk-3;

       trackL(1:k2-k1+1,2*ii-1:2*ii) = track(k1:k2,2*row-1:2*row);
   end

end


thrhd = 10/180*pi;
for ii = 1:size(trackL,2)./2
   jj = 1;
   sat = 0;
   temp = zeros(3,1);
   while ~sat 
       temp(:) = trackL(jj:jj+2,2*ii-1);
       if abs(temp(1)-temp(2))<=thrhd && abs(temp(2)-temp(3))<=thrhd
           sat = 1;
       else
           jj = jj +1;
       end      
   end
   
   rtrack(1,2*ii-1:2*ii) = mean(trackL(jj:jj+2,2*ii-1:2*ii));

   
   sat1 = 0;
   cnt = jj;
   while ~sat1
      temp1 = trackL(jj,2*ii-1);  
      if abs(temp1-rtrack(1,2*ii-1))<= thrhd && trackL(jj,2*ii)>0
          jj = jj+1;
      else
          sat1 = 1;
      end
   end
   rtrack(2,2*ii-1) = jj-cnt+1;
   rtrack(2,2*ii) = sum(trackL(cnt:jj,2*ii));
   
end
