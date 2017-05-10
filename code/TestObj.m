classdef TestObj < handle
    
   properties (SetAccess = private)
      Ed
      Ed0
      EdPts
      EdPts0
      Brch
      Brch0
      BrchPts
      BrchPts0
      fig
      BW
      Skel
      Pathway
      trajectory
   end
   
   methods
       % class constructor
       function this = TestObj()
           close all
           clc
           cd('C:\Users\lhuang28\Desktop\SimuVideo2 T\Alg1')
           open('figT.fig');
           this.fig = getframe;
           bw = imbinarize(this.fig.cdata(:,:,1),0.9);
           this.BW = imresize(bw,1); 
           this.Skel = bwmorph(this.BW,'skel',Inf);
           this.Brch = bwmorph(this.Skel,'branchpoints');
           [row col] = find(this.Brch);
           this.BrchPts = [row col]; 
           this.Ed = bwmorph(this.Skel,'endpoints');
           [row col] = find(this.Ed);
           this.EdPts = [row col];
           figure
           imshow(this.Skel);
           hold on
           %plot(BrchPts(:,2),BrchPts(:,1),'r.');
           f1 = scatter(this.BrchPts(:,2),this.BrchPts(:,1),20,'filled');
           f1.CData = [1 0 0];
           f2 = scatter(this.EdPts(:,2),this.EdPts(:,1),20,'filled');
           f2.CData = [0 0 1];
           %plot(EdPts(:,2),EdPts(:,1),'b.');
           this.Planning;
           this.preprocess;
           this.track;
       end
       
       function Planning(this)
            ArraySize = sum(this.Skel(:));
            this.Pathway = uint16(this.Skel);
            cost = 1e2;
            this.Pathway(245,273)=cost;
            frtr = [245,273];
            rg = 1;
            ii = 1;
            jj = 1;

            while ~isempty(frtr) %(jj<ArraySize)
                temp = this.Pathway(-1+frtr(ii,1):1+frtr(ii,1),-1+frtr(ii,2):1+frtr(ii,2));
                [row col] = find(temp==1);
                row = row-2;
                col = col-2;
                cnt = numel(row);
                if cnt>0
                    this.Pathway(-1+frtr(ii,1):1+frtr(ii,1),-1+frtr(ii,2):1+frtr(ii,2))=...
                        -this.Pathway(-1+frtr(ii,1):1+frtr(ii,1),-1+frtr(ii,2):1+frtr(ii,2)).*uint16(blkdiag(0,1,0))+...
                        uint16(temp ==1).*this.Pathway(frtr(ii,1),frtr(ii,2))...
                        +this.Pathway(-1+frtr(ii,1):1+frtr(ii,1),-1+frtr(ii,2):1+frtr(ii,2));

                    cost = cost +1;
                    for ij = 1:cnt
                        frtr = [frtr; row(ij)+frtr(ii,1), col(ij)+frtr(ii,2)];     
                    end
                end
                frtr(ii,:) = [];
                jj = jj + cnt;
            end

           
           
       end
       
       function preprocess(this)
          this.Ed0 = this.Ed;
          this.Brch0 = this.Brch;
          for ii = 1:length(this.EdPts)
                temp1 = this.EdPts(ii,1:2);
                temp = temp1;
                cost1 = this.Pathway(temp(1),temp(2));
                while(this.Brch(temp(1),temp(2))==0)
                    cost = this.Pathway(temp(1),temp(2));
                    nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                    nb(2,2)=0;
                    if ~isempty(find(nb==cost)) && ~isempty(find(nb==cost-1)) && ~isempty(find(nb==cost+1)) 
                        [row col] = find(nb==cost);
                        row = row-2;
                        col = col-2;
                        temp = [row(1)+temp(1),col(1)+temp(2)];
                        cost2 = this.Pathway(temp(1),temp(2));
                    else
                        [row col] = find(nb==cost-1);
                        row = row-2;
                        col = col-2;
                        temp = [row(1)+temp(1),col(1)+temp(2)];
                        cost2 = this.Pathway(temp(1),temp(2));
                    end
                end
               if (cost1-cost2<50)

                    this.Brch0(temp(1),temp(2))=0;
                    this.Ed0(temp(1),temp(2))=1;
                    this.Ed0(temp1(1),temp1(2))=0;
               end


            end

            [row col] = find(this.Brch0);
            this.BrchPts0 = [row col];
            [row col] = find(this.Ed0);
            this.EdPts0 = [row col];
           
       end
       
       function track(this)
          
            this.trajectory = [];
            for ii =1:length(this.EdPts0)
                temp1 = this.EdPts0(ii,1:2);
                temp2 = temp1;
                cost1 = this.Pathway(temp2(1),temp2(2));
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
                            while(cost_diff<26 && this.Brch0(temp2(1),temp2(2))==0)
                                nb = this.Pathway(-1+temp2(1):1+temp2(1),-1+temp2(2):1+temp2(2));
                                nb(2,2)=0;
                                if ~isempty(find(nb==cost)) 
                                    [row col] = find(nb==cost);
                                    row = row-2;
                                    col = col-2;
                                    temp = [row(1)+temp2(1),col(1)+temp2(2)];
                                    if this.Brch0(temp(1),temp(2))==0
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
                                cost2 = this.Pathway(temp2(1),temp2(2));
                                cost_diff = cost1-cost2;
                                cost = cost-1;
                            end
                        else
                             while cost_diff<26
                                nb = this.Pathway(-1+temp2(1):1+temp2(1),-1+temp2(2):1+temp2(2));
                                nb(2,2)=0;
                                if ~isempty(find(nb==cost)) 
                                    [row col] = find(nb==cost);
                                    row = row-2;
                                    col = col-2;
                                    temp = [row(1)+temp2(1),col(1)+temp2(2)];
                                    if this.Brch0(temp(1),temp(2))==0
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
                                cost2 = this.Pathway(temp2(1),temp2(2));
                                cost_diff = cost1-cost2;
                                cost = cost-1;
                            end

                        end

                        angle = atan2(double(temp2(1)-temp1(1)),double(temp2(2)-temp1(2)));
                        this.trajectory = [this.trajectory;angle, double(cost_diff)];
                    end
                    if this.Brch0(temp2(1),temp2(2))==1
                        flag = flag + 1/2;
                    else
                        if flag > 0
                            flag = flag + 1/2;
                        end
                    end
                end


                this.trajectory = [this.trajectory;NaN, NaN];

            end 
       end
   end
   
end