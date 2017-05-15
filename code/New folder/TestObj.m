classdef TestObj < handle
    
   properties (SetAccess = private)
      Ed
      Ed0
      EdPts
      EdPts0
      sdEdPts
      Brch
      Brch0
      BrchPts
      BrchPts0
      BEdist
      BBdist
      fig
      BW
      channel_width
      EBW
      Skel
      goal
      max_step
      ConnMat
      Pathway
      Path
      track
      mission_increment
      RegionID
      RegionVal
      trackL
      rtrack
      localSkel
      localPathway
      tmpPathway
      localPath
      tar
      tar_brch
      tar_brch_ID
      local_track
      distbt
      localmpath
      loc
      freespace
      localfsp
      medial
      basescore
      Expt
      row0
      col0
      rank0
   end
    
   methods
       %% class constructor
       function this = TestObj(map)
           
           %% Load Map
           testmap = strcat(map.name,'.fig');
           open(testmap);
           this.fig = getframe;
           this.channel_width = map.channel_width;
           this.mission_increment = map.mission_increment;
           bw = imbinarize(this.fig.cdata(:,:,1),0.9);      % RGB --> binary image
           this.BW = imresize(bw,map.magnification);        % Resize the image for smooth path
           this.max_step = map.max_step;                    % trajectory vector max magnitude
           % Augment the size of the binary image along boundaries
           [l10,l20] = size(this.BW);
           l1 = uint16(l10*1.2);
           l2 = uint16(l20*1.2);
           
           BW0 = uint16(zeros(l1,l2));
           BW0(uint16((l1-l10)/2):uint16((l1-l10)/2)-1+l10,uint16((l2-l20)/2):l20-1+uint16((l2-l20)/2))=this.BW;
           this.BW = logical(BW0);
           
           this.EBW = ~edge(this.BW,'canny');           
           this.Skel = bwmorph(this.BW,'skel',Inf);         % Obtain medial-axis skeleton of the binary image, targets are marked as 1, o.w. 0
           this.Brch = bwmorph(this.Skel,'branchpoints');   % Obtain branch map in the medial-axis map, targets are marked as 1, o.w. 0 
           [row, col] = find(this.Brch);                     
           this.BrchPts = [row, col];                        % Branch pts coordinates
           this.Ed = bwmorph(this.Skel,'endpoints');        % Obtain end-point map in the medial-axis map, targets are marked as 1, o.w. 0
           [row, col] = find(this.Ed);
           this.EdPts = [row, col];                          % End pts coordinates
           this.goal = map.goal_loc;
           this.BEdist = map.distance_threshold(1);
           this.BBdist = map.distance_threshold(2);
           %% Visualize the map
% % %            figure
% % %            imshow(this.Skel);
% % %            hold on
% % %            f1 = scatter(this.BrchPts(:,2),this.BrchPts(:,1),20,'filled');
% % %            f1.CData = [1 0 0];
% % %            f2 = scatter(this.EdPts(:,2),this.EdPts(:,1),20,'filled');
% % %            f2.CData = [0 0 1];
% % %            
           
           %% Part I: Map Preprocessing
           MapProcess1(this);
           MapProcess2(this);
           MapProcess3(this);
           MapProcess4(this);
           %% Part II: Tracjetory Vectors Generation
           MapProcess5(this);
           MapProcess6(this);
           %% Part III: Local Map Evaluation
         
           this.rank0 = 0;
           
           for ii = 6:numel(this.RegionID)
               
               this.rank0 = ii;
               LocalMapProc1(this);
               [row, col] = find(mod(this.tmpPathway,round(this.channel_width/3))==0);
               for jj = 1:numel(row)
                   this.col0 = jj;
                   
                   for kk = 1:10
                        this.row0 = kk+2; 
                        LocalMapProc1(this,jj,ii);
                   end
               end
           end
%            LocalMapProc2(this);
%            LocalMapProc3(this);
%            GlobalControl(this);
%            LocalControl(this); 
       end
       
       %% Assign cost (from the goal)to nodes
       function MapProcess1(this)   
            this.Pathway = uint16(this.Skel);                   % 'Pathway' is a cost-to-go value map
            temp = this.Skel == 1;
            
            this.Path = zeros(sum(temp(:)),6);
            NodeCnt = 1;
            parentID = 0;

            this.Pathway(this.goal(1),this.goal(2))=1e2;       % goal cost
            frtr = this.goal;                                   % BFS frontier
            this.sdEdPts = [];
            this.Path(NodeCnt,1:2) = this.goal;
            this.Path(NodeCnt,4) = 1e2;
            
            while ~isempty(frtr) 
                % Extract a 3x3 neighbor of the current frontier node from 'Pathway'  
                nb = this.Pathway(-1+frtr(1,1):1+frtr(1,1),-1+frtr(1,2):1+frtr(1,2));     
                [row, col] = find(nb==1);      % Search for unexplored neighbor pts that belongs to medial-axie skeleton map          
                row = row-2;
                col = col-2;
                
                if numel(row)>0
                    % Assign cost value to the new frontier
                    this.Pathway(-1+frtr(1,1):1+frtr(1,1),-1+frtr(1,2):1+frtr(1,2))=...
                        -this.Pathway(-1+frtr(1,1):1+frtr(1,1),-1+frtr(1,2):1+frtr(1,2)).*uint16(blkdiag(0,1,0))+...
                        uint16(nb ==1).*this.Pathway(frtr(1,1),frtr(1,2))...
                        +this.Pathway(-1+frtr(1,1):1+frtr(1,1),-1+frtr(1,2):1+frtr(1,2));
                    parentID = find(this.Path(1:NodeCnt,1)==frtr(1,1) & this.Path(1:NodeCnt,2)==frtr(1,2));
                    % Add new frontier nodes to the var 'frtr'
                    for ij = 1:numel(row)       
                        frtr = [frtr; row(ij)+frtr(1,1), col(ij)+frtr(1,2)];
                        NodeCnt = NodeCnt +1;
                        this.Path(NodeCnt,1:2) = frtr(end,1:2);
                        this.Path(NodeCnt,3) = parentID;
                        this.Path(NodeCnt,4) = this.Pathway(frtr(1,1),frtr(1,2))+1;
                    end
                else
                    % sdEdPts stands for pseudo end pts, which is not an
                    % actual end pts, but a pts where you can get to the
                    % goal with more than one path at the same cost
                    if this.Ed(frtr(1,1),frtr(1,2))==0
                        this.sdEdPts = [this.sdEdPts;frtr(1,1:2)]; 
                    end
                end
                % The current node is no more at the frontier
                frtr(1,:) = []; 
            end
            
            
            if ~isempty(this.sdEdPts)
                this.sdEdPts = this.sdEdPts(find(diag(this.Brch(this.sdEdPts(:,1),this.sdEdPts(:,2)))==0),:);
            end
       end
          
       %% Eliminate Fake End Pts
       function MapProcess2(this)
             % Make sure the pts with the largest cost are end pts 
             [row, col] = find(this.Pathway==max(this.Pathway(:)));
              for jj = 1:numel(row)
                    if this.Ed(row(jj),col(jj))==0
                       this.Ed(row(jj),col(jj))=1;
                    end
              end
              
              % Pass teh info of branch pts and end pts to process vars
              % e. g. this.Ed is the original data
              % this.Ed0 is the data for processing
              
              this.Ed0 = this.Ed;
              this.Brch0 = this.Brch;
              
              % Make sure the goal is an end pts and a branch pts
              this.Ed0(this.goal(1),this.goal(2)) = 1;
              this.Brch0(this.goal(1),this.goal(2)) = 1;

              [row, col] = find(this.Ed0==1);
              this.EdPts0 = [row, col];

              [row, col] = find(this.Brch0==1);
              this.BrchPts0 = [row, col];
                
              % This for-loop says that if the distance (cost) of a branch
              % pts and the nearest end pts is less than the threshold (this.BEdist),
              % then this branch pts is regarded as an end pts instead of a
              % branch pts.
              for ii = 1:length(this.EdPts)
                    temp1 = this.EdPts(ii,1:2);                 % Pick an end point
                    temp = temp1;                               % Intermediate pt for BFS search
                    cost1 = this.Pathway(temp1(1),temp1(2));    % This end point cost
                    cost = cost1;
                    tempRcd = temp;
                    while (cost1-cost < this.BEdist) && cost > 100
                        
                        % Extract the neighbors around the intermediate pt
                        % from the cost map 'Pathway'
                        nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb(2,2)=0;
                        
                        % Extract the neighbors around the intermediate pt
                        % from the Branch map 'Brch'
                        nb0 = this.Brch(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb0(2,2)=0;
                        
                        % Search for branch pts
                        if ~isempty(find(nb0==1,1))
                            [row, col] = find(nb0==1);
                            row = row-2;
                            col = col-2;
                            if this.Pathway(row(1)+temp(1),col(1)+temp(2)) > cost
                               [row, col] = find(nb==cost-1);
                               row = row-2;
                               col = col-2;
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                               tempRcd = [tempRcd;temp];
                            else
                               % The branch pt is within the threshold distance, and it
                               % will be regarded as an end pt.    
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                               this.Brch0(temp(1),temp(2))=0;
                               this.Ed0(temp(1),temp(2))=1;
                               this.Ed0(temp1(1),temp1(2))=0;
                               for i_tempRcd = 1:size(tempRcd,1)
                                    this.Skel(tempRcd(i_tempRcd,1),tempRcd(i_tempRcd,2)) = 0;
                                    this.Pathway(tempRcd(i_tempRcd,1),tempRcd(i_tempRcd,2)) = 0;
                               end
                                Idx = knnsearch(this.Path(:,1:2),tempRcd);
                                this.Path(Idx,:) = NaN;
                               break;
                            end
                        else
                            [row, col] = find(nb==cost-1);
                            row = row-2;
                            col = col-2;
                            temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                            tempRcd = [tempRcd;temp];
                        end
                        
                        cost = this.Pathway(temp(1),temp(2));       % Intermediate pt cost
                        
                   end
                    
                   
              end
       
              % Update end pts list
              [row, col] = find(this.Ed0);
              this.EdPts0 = [row, col];
           
              % This for-loop eliminates any possible 'end pts' b/w one end
              % pt and the goal
              for ii = 1:size(this.EdPts0,1)
                    temp = this.EdPts0(ii,:);               % Obtain coordinates of this end pt                   
                    
                    if this.Ed0(temp(1),temp(2))==0         % Is this an end pt? We are updating end pts info
                        continue
                    end

                    cost = this.Pathway(temp(1),temp(2));
                    while cost > 100    
                        % Extract the neighbors around the intermediate pt
                        % from the cost map 'Pathway'
                        nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb(2,2)=0;
                        
                        % Extract the neighbors around the intermediate pt
                        % from the Branch map 'Brch'
                        nb0 = this.Ed0(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb0(2,2)=0;
                        
                                                % Search for branch pts
                        if ~isempty(find(nb0==1,1))
                            [row, col] = find(nb0==1);
                            row = row-2;
                            col = col-2;
                            if this.Pathway(row(1)+temp(1),col(1)+temp(2)) > cost
                               [row, col] = find(nb==cost-1);
                               row = row-2;
                               col = col-2;
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt

                            else
                               % This 'end pt' is b/w an end pt and the goal, and it
                               % will NOT be regarded as an end pt.    
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                               this.Ed0(temp(1),temp(2))=0;
                            end
                        else
                            [row, col] = find(nb==cost-1);
                            row = row-2;
                            col = col-2;
                            temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                        end
                        
                        cost = this.Pathway(temp(1),temp(2));       % Intermediate pt cost

                    end 
              end

            
            
       end
        
       %% Eliminate Fake Branch Pts
       function MapProcess3(this)
                
              % This for-loop says that if the distance (cost) of two nearby branch
              % pts is less than the threshold (this.BBdist),
              % then the farther branch pt will NOT be regarded as a branch pts
              for ii = 1:size(this.BrchPts0,1)
                    temp1 = this.BrchPts0(ii,1:2);                 % Pick an end point
                    temp = temp1;                               % Intermediate pt for BFS search
                    cost1 = this.Pathway(temp1(1),temp1(2));    % This end point cost
                    cost = cost1;
                    
                    while (cost1-cost < this.BBdist) && cost >100
                        
                        % Extract the neighbors around the intermediate pt
                        % from the cost map 'Pathway'
                        nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb(2,2)=0;
                        
                        % Extract the neighbors around the intermediate pt
                        % from the Branch map 'Brch'
                        nb0 = this.Brch0(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb0(2,2)=0;
                        
                        % Search for branch pts
                        if ~isempty(find(nb0==1,1))
                            [row, col] = find(nb0==1);
                            row = row-2;
                            col = col-2;
                            if this.Pathway(row(1)+temp(1),col(1)+temp(2)) > cost
                               [row, col] = find(nb==cost-1);
                               row = row-2;
                               col = col-2;
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt

                            else
                               % The branch pt is within the threshold distance, and it
                               % will be regarded as an end pt.    
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                               this.Brch0(temp1(1),temp1(2))=0;

                            end
                        else
                            [row, col] = find(nb==cost-1);
                            row = row-2;
                            col = col-2;
                            temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                        end
                        
                        cost = this.Pathway(temp(1),temp(2));       % Intermediate pt cost
                        
                   end
                    
                   
              end    
              
            this.Brch0(this.goal(1),this.goal(2))=1;
            [row, col] = find(this.Brch0);
            this.BrchPts0 = [row, col];
            this.Ed0(this.goal(1),this.goal(2)) = 1;
            [row, col] = find(this.Ed0);
            this.EdPts0 = [row, col];
            
% % %             figure
% % %             imshow(~(~this.BW+this.Skel)); hold on
% % %             
% % %             
% % %             %plot(BrchPts(:,2),BrchPts(:,1),'r.');
% % %             f1 = scatter(this.BrchPts0(:,2),this.BrchPts0(:,1),20,'filled');
% % %             f1.CData = [1 0 0];
% % %             f2 = scatter(this.EdPts0(:,2),this.EdPts0(:,1),20,'filled');
% % %             f2.CData = [0 0 1];
% % %             if ~isempty(this.sdEdPts)
% % %                 f3 = scatter(this.sdEdPts(:,2),this.sdEdPts(:,1),5,'filled');
% % %                 f3.CData = [0,0,0.8];
% % %             end

       end
            
       %% Pair junction (branch pts) with end pts
       function MapProcess4(this)
           % Note this.track will change row size later.
           % The j-th col indicates trajectories (in the form of vectors) lead to
           % the j-th branch pt (sorted by cost-to-go ('descend')).
           % In one col, trajectories start from different end pt/ branch
           % pt are sepqrated by "NaN". There are multiple trajectories
           % lead to one branch pt
           
           this.track = zeros(1,2*size(this.BrchPts0,1));
           % Note this.trackL will change row/col size later.
           % The j-th col indicates a trajectory (by tracjectory vectors)
           % from j-th branch pt (sorted by cost-to-go ('descend')) to the next branch pt (junction)
           % Compared to a col in this.track, a col in this.trackL only contains 1 traj. 
           this.trackL = zeros(1,2*size(this.BrchPts0,1));
           
           this.rtrack = this.trackL;
           
            % Sort the region by cost-to-go from the branch pts
           [this.RegionVal, this.RegionID] = sort(diag(this.Pathway(this.BrchPts0(:,1),this.BrchPts0(:,2))),'Descend');
           
           this.ConnMat = zeros(size(this.BrchPts0,1),11);
           
          % Pair branch pts w/ branch pts        
            for ii = 1:numel(this.RegionID)
                temp1 = this.BrchPts0(this.RegionID(ii),1:2);    % Pick an end point
                temp = temp1;                               % Intermediate pt for BFS search
                cost1 = this.RegionVal(ii);                      % This end point cost
                cost = cost1;

                while cost >100

                    % Extract the neighbors around the intermediate pt
                    % from the cost map 'Pathway'
                    nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                    nb(2,2)=0;

                    % Extract the neighbors around the intermediate pt
                    % from the Branch map 'Brch'
                    nb0 = this.Brch0(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                    nb0(2,2)=0;

                    % Search for branch pts
                    if ~isempty(find(nb0==1,1))
                        [row, col] = find(nb0==1);
                        row = row-2;
                        col = col-2;
                        % If the branch pt has higher cost than the intermediate pt "temp", skip it
                        if this.Pathway(row(1)+temp(1),col(1)+temp(2)) > cost
                           [row, col] = find(nb==cost-1);
                           row = row-2;
                           col = col-2;
                           temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt

                        else
                           % If we reach a branch pt, break the while loop
                           % Assign the starting pt to this branch pt.
                           temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                           ID = find(this.BrchPts0(:,1)==temp(1) & this.BrchPts0(:,2)==temp(2));
                           this.ConnMat(ID,1) = this.ConnMat(ID,1)+1;
                           this.ConnMat(ID,2*this.ConnMat(ID,1):2*this.ConnMat(ID,1)+1) = temp1;
                           cost = this.Pathway(temp(1),temp(2));
                           break;
                        end
                    else
                        % No branch pts found, keep looking for
                        % intermediate pt of less cost: (cost - 1)
                        [row, col] = find(nb==cost-1);
                        row = row-2;
                        col = col-2;
                        temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                    end

                    cost = this.Pathway(temp(1),temp(2));       % Intermediate pt cost

                end


            end
            
            % Sort end pts by cost-to-go
            [EndVal, EndID] = sort(diag(this.Pathway(this.EdPts0(:,1),this.EdPts0(:,2))),'Descend');
            
            % Pair end pts w/ branch pts        
            for ii = 1:numel(EndID)
                temp1 = this.EdPts0(EndID(ii),1:2);         % Pick an end point
                temp = temp1;                               % Intermediate pt for BFS search
                cost1 = EndVal(ii);                         % This end point cost
                cost = cost1;

                while cost > 100

                    % Extract the neighbors around the intermediate pt
                    % from the cost map 'Pathway'
                    nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                    nb(2,2)=0;

                    % Extract the neighbors around the intermediate pt
                    % from the Branch map 'Brch'
                    nb0 = this.Brch0(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                    nb0(2,2)=0;

                    % Search for branch pts
                    if ~isempty(find(nb0==1,1))
                        [row, col] = find(nb0==1);
                        row = row-2;
                        col = col-2;
                        
                        % If the branch pt has higher cost than the intermediate pt "temp", skip it
                        if this.Pathway(row(1)+temp(1),col(1)+temp(2)) > cost
                           [row, col] = find(nb==cost-1);
                           row = row-2;
                           col = col-2;
                           temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt

                        else
                           % If we reach a branch pt, break the while loop
                           % Assign the starting pt to this branch pt. 
                           temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                           ID = find(this.BrchPts0(:,1)==temp(1) & this.BrchPts0(:,2)==temp(2));
                           this.ConnMat(ID,1) = this.ConnMat(ID,1)+1;
                           this.ConnMat(ID,2*this.ConnMat(ID,1):2*this.ConnMat(ID,1)+1) = temp1;
                           cost = this.Pathway(temp(1),temp(2));
                           break;
                        end
                    else
                        % No branch pts found, keep looking for
                        % intermediate pt of less cost: (cost - 1)
                        [row, col] = find(nb==cost-1);
                        row = row-2;
                        col = col-2;
                        temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                    end
                    cost = this.Pathway(temp(1),temp(2));       % Intermediate pt cost

                end


            end
            
            % If there exists pseudo end pts (the pts cannot be found by bwmorph)
            if ~isempty(this.sdEdPts)
                % Sort sdEdPts
                [EndVal, EndID] = sort(diag(this.Pathway(this.sdEdPts(:,1),this.sdEdPts(:,2))),'Descend');
            
                % Pair pseudo end pts w/ branch pts        
                for ii = 1:numel(EndID)
                    temp1 = this.sdEdPts(EndID(ii),1:2);                 % Pick an end point
                    temp = temp1;                               % Intermediate pt for BFS search
                    cost1 = EndVal(ii);    % This end point cost
                    cost = cost1;

                    while cost > 100

                        % Extract the neighbors around the intermediate pt
                        % from the cost map 'Pathway'
                        nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb(2,2)=0;

                        % Extract the neighbors around the intermediate pt
                        % from the Branch map 'Brch'
                        nb0 = this.Brch0(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb0(2,2)=0;

                        % Search for branch pts
                        if ~isempty(find(nb0==1,1))
                            [row, col] = find(nb0==1);
                            row = row-2;
                            col = col-2;
                            % If the branch pt has higher cost than the intermediate pt "temp", skip it
                            if this.Pathway(row(1)+temp(1),col(1)+temp(2)) > cost
                               [row, col] = find(nb==cost-1);
                               row = row-2;
                               col = col-2;
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt

                            else
                               % If we reach a branch pt, break the while loop
                               % Assign the starting pt to this branch pt. 
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                               ID = find(this.BrchPts0(:,1)==temp(1) & this.BrchPts0(:,2)==temp(2));
                               this.ConnMat(ID,1) = this.ConnMat(ID,1)+1;
                               this.ConnMat(ID,2*this.ConnMat(ID,1):2*this.ConnMat(ID,1)+1) = temp1;
                               cost = this.Pathway(temp(1),temp(2));
                               break;
                            end
                        else
                        % No branch pts found, keep looking for
                        % intermediate pt of less cost: (cost - 1)
                            [row, col] = find(nb==cost-1);
                            row = row-2;
                            col = col-2;
                            temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                        end

                        cost = this.Pathway(temp(1),temp(2));       % Intermediate pt cost

                    end


                end
            
            end
            
 

       end
       
       %% Trajectory Vectors Generation
       function MapProcess5(this)
            temp1 = this.ConnMat(1,2*1:2*1+1);
            % Initialization of target-motion Demo
%             h = scatter(temp1(2),temp1(1),80,'filled');
%             h.CData = [1 0 0];

            % Loop all the branch pts (junctions)
            for jj = 1:size(this.ConnMat,1)
                kk = 1;
                % Loop all the end pts/branch pts connecting to the current
                % branch pt (junction)
                for ii =1:this.ConnMat(this.RegionID(jj),1)
                    % Start from the selected end/branch pt and head the
                    % junction.
                    temp1 = this.ConnMat(this.RegionID(jj),2*ii:2*ii+1);
                    ks = [];
                    
                    % If this trail connects two branch pts, we record the
                    % starting pt and save into this.trackL latter
                    if this.Brch0(temp1(1),temp1(2)) == 1 && this.Ed0(temp1(1),temp1(2)) ==0
                        ks = kk;
                    end
                    
                    
                    temp = temp1;
                    cost1 = this.Pathway(temp1(1),temp1(2));
                    cost = cost1;
                    
                    mission_flag = 0;
                    % mission_flag is an indicator for the trail
                    % exploration. Starting from a end pt/branch pt towards
                    % the junction, if we reach the junction, mission_flag
                    % = mission_flag + mission_increment (now
                    % mission_flag<1), and we keep going for (1/mission_increment-1)*max_step 
                    % before we stop, because we need to send the robot
                    % into the next region instead of stopping at the
                    % junction
                    while mission_flag<=1 && cost>100
                        temp1 = temp;
                        cost1 = cost;                         
                        % While the trajectory vector (from temp1 to temp)
                        % is less this.max_step
                        while cost1-cost < this.max_step
                            % Extract the 3*3 neighborhood of the current location "temp" in "Pathway"
                            nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                            nb(2,2)=0;
                            cost = this.Pathway(temp(1),temp(2));
                            % Extract the 3*3 neighborhood of the current location "temp" in "Brch0"
                            nb0 = this.Brch0(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                            nb0(2,2)=0;

                            % If there's a branch pt (junction) around
                            if ~isempty(find(nb0==1,1)) 
                                [row, col] = find(nb0==1);
                                row = row-2;
                                col = col-2;
                                % And this branch pt is behind the
                                % current location, we skip it and use
                                % gradient descent to look for the next
                                % node
                                if this.Pathway(row(1)+temp(1),col(1)+temp(2)) > cost
                                   [row, col] = find(nb==cost-1);
                                   row = row-2;
                                   col = col-2;
                                   temp = [row(1)+temp(1),col(1)+temp(2)];
                                else
                                % else means that we have reached the first junction, 
                                % so stop increasing the trajectory vector length 
                                % and update mission_flag --> mission_flag >0                                
                                   temp = [row(1)+temp(1),col(1)+temp(2)];
                                   mission_flag = mission_flag + this.mission_increment;
                                   cost = this.Pathway(temp(1),temp(2));
                                   break;
                                end
                            else
                                 % else we use gradient descent to look for
                                % the next node
                                [row, col] = find(nb==cost-1);
                                row = row-2;
                                col = col-2;
                                temp = [row(1)+temp(1),col(1)+temp(2)];
                            end
                            cost = this.Pathway(temp(1),temp(2));
                            % Update the target-motion Demo
%                             h.XData = temp(2);
%                             h.YData = temp(1);
%                             pause(0.001);
                        end
                        % If we have reached a junction before,
                        % mission_flag>0, we keep going for
                        % (1./mission_increment-1)*max_step before we stop to look
                        % at another end pt.
                        mission_flag = mission_flag + ceil(mission_flag)*this.mission_increment;
                        % calculate the trajectory vector orientation
                        angle = atan2(double(temp(1)-temp1(1)),double(temp(2)-temp1(2)));

                        if cost1 - cost > 0
                            % Update trajectory vectors with its orientation and length 
                            this.track(kk,2*jj-1:2*jj) = [angle, double(cost1-cost)];
                            kk = kk+1;
                        end
                    end
                  
                    % From ks to kk we record a trajectory from a branch pt
                    % to another branch pt (the other situation is from an end pt to a branch pt)
                    % 
                    if ~isempty(ks)
                       % find the ID of this starting branch pt
                       ID = find(this.BrchPts0(:,1)==this.ConnMat(this.RegionID(jj),2*ii) ...
                        & this.BrchPts0(:,2)==this.ConnMat(this.RegionID(jj),2*ii+1));
                      % record this trajectory in this.trackL for latter
                      % use
                        this.trackL(1:kk-1/this.mission_increment-ks+1, ...
                           2*this.RegionID(ID)-1:2*this.RegionID(ID))= ...
                           this.track(ks:kk-1/this.mission_increment,2*jj-1:2*jj);
                        ks = [];
                    end
                    % After we reach the 1st branch pt after we start, we
                    % step a little bit further then stop. Then we record
                    % this trajectory in this.track, and use "NaN" to separate it from 
                    % the next trajectory leads to the same branch pt
                    this.track(kk,2*jj-1:2*jj) = [NaN, NaN];
                    kk = kk+1;

                end

            end
            
            this.trackL(end+1,:) = 0;
       end
      
       %% Trajectory Stats    
       function MapProcess6(this)
           % Trajectory smoothness threshold
            thrhd = 10/180*pi;
            for ii = 1:size(this.BrchPts0,1)
               jj = 1;
               flag = 0;
               temp = zeros(3,1);
               while ~flag 
                   temp(:) = this.trackL(jj:jj+2,2*ii-1);
                   if abs(temp(1)-temp(2))<=thrhd && abs(temp(2)-temp(3))<=thrhd
                       flag = 1;
                   else
                       jj = jj +1;
                   end      
               end

               this.rtrack(1,2*ii-1:2*ii) = mean(this.trackL(jj:jj+2,2*ii-1:2*ii));
               flag1 = 0;
               init = jj;
               while ~flag1
                  temp = this.trackL(jj,2*ii-1);  
                  if abs(temp-this.rtrack(1,2*ii-1))<= thrhd && this.trackL(jj,2*ii)>0
                      jj = jj+1;
                  else
                      flag1 = 1;
                  end
               end
               this.rtrack(2,2*ii-1) = jj-init+1;
               this.rtrack(2,2*ii) = sum(this.trackL(init:jj,2*ii));

            end
       
       end
       
       
       function LocalMapProc1(this,arg1,arg2)
         
           % Pick a random reion
           this.tar_brch = [];
         
           if nargin < 2
                this.tar_brch_ID = this.rank0;
                this.tar_brch = this.BrchPts0(this.RegionID(this.tar_brch_ID),1:2);
           else
               this.tar_brch_ID = arg2;
               this.tar_brch = this.BrchPts0(this.RegionID(this.tar_brch_ID),1:2);
           end
           
           this.localSkel = this.Skel.*0;
           this.localSkel(this.tar_brch(1),this.tar_brch(2)) = 1;
           
           % Loop all end pts of the targeted region tar_brch 
           % plus the branch pt (junction node)
           for ii = 1:this.ConnMat(this.RegionID(this.tar_brch_ID),1)+1  
                    if ii == this.ConnMat(this.RegionID(this.tar_brch_ID),1)+1
                        % select the branch pt (the junction)
                        temp = this.BrchPts0(this.RegionID(this.tar_brch_ID),1:2);
                    else
                        % select an end pt
                        temp = this.ConnMat(this.RegionID(this.tar_brch_ID),2*ii:2*ii+1);     
                    end
                    cost = this.Pathway(temp(1),temp(2));   % This end point cost
                    this.localSkel(temp(1),temp(2)) = 1;
                    
                    while cost >100
                        
                        % Extract the neighbors around the intermediate pt
                        % from the cost map 'Pathway'
                        nb = this.Pathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb(2,2)=0;
                        
                        % Extract the neighbors around the intermediate pt
                        % from the Branch map 'Brch'
                        nb0 = this.Brch0(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        nb0(2,2)=0;
                        
                        % Search for branch pts
                        if ~isempty(find(nb0==1,1))
                            [row, col] = find(nb0==1);
                            row = row-2;
                            col = col-2;
                            if this.Pathway(row(1)+temp(1),col(1)+temp(2)) > cost
                               [row, col] = find(nb==cost-1);
                               row = row-2;
                               col = col-2;
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt

                            else
                               % The branch pt is within the threshold distance, and it
                               % will be regarded as an end pt.    
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                               this.localSkel(temp(1),temp(2))=1;
                               break;

                            end
                        else
                            [row, col] = find(nb==cost);
                            if ~isempty(row)
                                row = row-2;
                                col = col-2;
                                this.localSkel(row(1)+temp(1),col(1)+temp(2))=1;
                            end
                            [row, col] = find(nb==cost-1);
                            row = row-2;
                            col = col-2;
                            temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                        end
                        this.localSkel(temp(1),temp(2))=1;
                        cost = this.Pathway(temp(1),temp(2));       % Intermediate pt cost
                        
                   end               
           
           end
           
           [row, col] = find(this.localSkel==1);
           localPts = [row, col];
           if nargin < 2 
                this.tar = this.ConnMat(this.RegionID(this.tar_brch_ID),2:3);
                this.basescore = [];
           else
            
               [row, col] = find(mod(this.tmpPathway,round(this.channel_width/4))==0);
               this.tar = [row(arg1),col(arg1)];
           end
           cost = 1e2;
           this.localPathway = uint16(this.localSkel);
           temp = this.localSkel==1;
           this.localPath = zeros(sum(temp(:)),6);
           this.localPathway(this.tar(1),this.tar(2)) = cost;
           NodeCnt = 1;
           this.localPath(NodeCnt,1:2) = this.tar;
           this.localPath(NodeCnt,4) = 1e2;
           % Frontier Node
           frtr = this.tar;
           while ~isempty(frtr)
                 
                    % Extract a 3x3 neighbor of the current frontier node from 'Pathway'  
                    nb = this.localPathway(-1+frtr(1,1):1+frtr(1,1),-1+frtr(1,2):1+frtr(1,2));     
                    [row, col] = find(nb==1);      % Search for unexplored neighbor pts that belongs to medial-axie skeleton map          
                    row = row-2;
                    col = col-2;

                    if numel(row)>0
                        % Assign cost value to the new frontier
                        this.localPathway(-1+frtr(1,1):1+frtr(1,1),-1+frtr(1,2):1+frtr(1,2))=...
                            -this.localPathway(-1+frtr(1,1):1+frtr(1,1),-1+frtr(1,2):1+frtr(1,2)).*uint16(blkdiag(0,1,0))+...
                            uint16(nb ==1).*this.localPathway(frtr(1,1),frtr(1,2))...
                            +this.localPathway(-1+frtr(1,1):1+frtr(1,1),-1+frtr(1,2):1+frtr(1,2));

                        parentID = find(this.localPath(1:NodeCnt,1)==frtr(1,1) & this.localPath(1:NodeCnt,2)==frtr(1,2));
                        % Add new frontier nodes to the var 'frtr'
                        for ij = 1:numel(row)       
                            frtr = [frtr; row(ij)+frtr(1,1), col(ij)+frtr(1,2)];
                            NodeCnt = NodeCnt + 1;
                           
                            this.localPath(NodeCnt,1:2) = frtr(end,1:2);
                            this.localPath(NodeCnt,3) = parentID;
                            this.localPath(NodeCnt,4) = this.localPathway(frtr(1,1),frtr(1,2))+1;
                        end
                    end
                    % The current node is no more at the frontier
                    frtr(1,:) = []; 
                
               
               
           end
           
           if nargin < 2
               this.tmpPathway = this.localPathway+1;
           end
               GlobalControl(this);
               LocalControl(this);  


           
       end
       
       function LocalMapProc2(this)
            this.distbt = zeros(this.ConnMat(this.RegionID(this.tar_brch_ID),1),5);
           
            figure             
            
            imshow(~(this.localSkel+~this.BW));
            hold on
            scatter(this.tar(2),this.tar(1),100,'filled','p');
            temp1 = this.ConnMat(this.RegionID(this.tar_brch_ID),2:3);
            temp = temp1;
            % Initialization of target-motion Demo
            hlc = scatter(temp(2),temp(1),35,'filled');
            hlc.CData = [0 1 0];

            % Loop all the end pts/branch pts connecting to the current
            % branch pt (junction)
            for ii =1:this.ConnMat(this.RegionID(this.tar_brch_ID),1)
                % Start from the selected end/branch pt and head the
                % junction.
                kk = 1;
                temp1 = this.ConnMat(this.RegionID(this.tar_brch_ID),2*ii:2*ii+1);
                this.distbt(ii,3) = this.localPathway(temp1(1),temp1(2))-100;
                this.distbt(ii,4) = this.Pathway(temp1(1),temp1(2))-this.Pathway(this.tar(1),this.tar(2));
                this.distbt(ii,1:2) = temp1;
                temp = temp1;
                cost1 = this.localPathway(temp1(1),temp1(2));
                cost = cost1;

                while cost > 100
                    temp1 = temp;
                    cost1 = cost;                         
                    % While the trajectory vector (from temp1 to temp)
                    % is less this.max_step
                    while cost1-cost < this.max_step
                        % Extract the 3*3 neighborhood of the current location "temp" in "Pathway"
                        nb = this.localPathway(-1+temp(1):1+temp(1),-1+temp(2):1+temp(2));
                        cost = this.localPathway(temp(1),temp(2));

                        % If there's a branch pt (junction) around
                        if ~isempty(find(cost==100,1)) 
                            break;
                        else
                            % else we use gradient descent to look for
                            % the next node
                            [row, col] = find(nb==cost-1);
                            row = row-2;
                            col = col-2;
                            temp = [row(1)+temp(1),col(1)+temp(2)];
                        end
                        cost = this.localPathway(temp(1),temp(2));

                    end
                    hlc.XData = temp(2);
                    hlc.YData = temp(1);
                    pause(0.0001)
                    % calculate the trajectory vector orientation
                    angle = atan2(double(temp(1)-temp1(1)),double(temp(2)-temp1(2)));

                    if cost1 - cost > 0
                        % Update trajectory vectors with its orientation and length 
                        this.local_track(kk,2*ii-1:2*ii) = [angle, double(cost1-cost)];
                        kk = kk+1;
                    end
                end
            end

            this.distbt(:,5) = this.distbt(:,3)==this.distbt(:,4);
            this.local_track(end+1,:) = 0;
       end
       
       function LocalMapProc3(this)
          this.localmpath = zeros(2,2);
          % Find the branch the target loc is in 
          this.localmpath(1,1) = find(this.distbt(:,5)==1);
          % problem: tar at different loc.
          
          
          
          % Check the branch depth from the end pt to the tar
          this.localmpath(1,2) = find(this.local_track(:,2*this.localmpath(1,1))==0)-1;
          if this.localmpath(1,2) <=3 
                this.localmpath(2,1) = mean(this.local_track(1:this.localmpath(1,2),2*this.localmpath(1,1)-1));
                this.localmpath(2,2) = sum(this.local_track(1:this.localmpath(1,2),2*this.localmpath(1,1)));
          else
                this.localmpath(2,1) = mean(this.local_track(this.localmpath(1,2)-3:this.localmpath(1,2),2*this.localmpath(1,1)-1));
                this.localmpath(2,2) = sum(this.local_track(this.localmpath(1,2)-3:this.localmpath(1,2),2*this.localmpath(1,1)));
                jj = this.localmpath(1,2)-1;
                thrhd = 10/180*pi;
                flag = 0;
                while ~flag
                   temp = this.local_track(jj,2*this.localmpath(1,1)-1:2*this.localmpath(1,1));
                   if abs(temp-this.localmpath(2,1))<thrhd
                       this.localmpath(2,2) = this.localmpath(2,2)+temp(2);
                       jj = jj - 1;
                       if jj < 1
                           flag = 1;
                       end
                   else
                       flag = 1;
                   end
                    
                end
          end
          
          CMreact = [];
          react = [];
          seq = perms(linspace(1,this.ConnMat(this.tar_brch_ID,1),this.ConnMat(this.tar_brch_ID,1)));
          for ii = 1:factorial(this.ConnMat(this.tar_brch_ID,1))
                mm = 2;
                for jj = 1:this.ConnMat(this.tar_brch_ID,1)
                    
                   kk = 1;
                   ll = seq(ii,jj);
                   while this.local_track(kk,2*ll) > 0
                       temp = YC(ii,mm-1)+cur_track(kk,2*ll).*cos(cur_track(kk,2*ll-1)-tar(1));
                       if temp >tar(3)
                            temp = tar(3);
                       end
                       Y(ii,mm) = temp-YC(ii,mm-1);
                       YC(ii,mm) = temp;
                       mm = mm + 1;
                       kk = kk + 1;
                   end
                end
          end



           
       end
       
       function GlobalControl(this)             
            temp = this.BrchPts0(this.RegionID,1:2);
            for ii = 1:length(temp)
               Idx = find(this.Path(:,1)==temp(ii,1) & this.Path(:,2)==temp(ii,2));
               this.Path(Idx,5) = 1;
               this.Path(Idx,6) = ii;
            end
            

            temp = find(this.Path(:,5)==1);
            Idx = rangesearch(this.Path(:,1:2),this.Path(temp,1:2),2);
            for ii = 1:size(temp,1)
                this.Path(Idx{ii},6) = this.Path(temp(ii),6);
                this.Path(Idx{ii},5) = 0.5;
            end
            
            for ii = 1:length(this.Path)
               if this.Path(ii,6) ==0
                   this.Path(ii,6) = this.Path(this.Path(ii,3),6);
               end
            end
            
            [row, col] =  find(this.BW==1);
            this.freespace = [row, col];
            this.medial = this.Path(:,1:2);
            fspIdx = knnsearch(this.medial,this.freespace);
            this.freespace(:,4) = fspIdx;
            for ii = 1:numel(this.RegionID)
              temp = [];  
              temp1 = find(this.Path(:,6)==ii);
              for jj = 1:length(temp1)
                    temp2 = find(fspIdx == temp1(jj));
                    temp = [temp;temp2];
              end
              this.freespace(temp,3) = ii;
            end
            
% % %             figure
% % %             imshow(this.BW);
% % %             hold on
% % %             for ii = 1:numel(this.RegionID)
% % %                 temp = find(this.freespace(:,3)==ii);
% % %                 scatter(this.freespace(temp,2),this.freespace(temp,1),10,'filled');
% % %             end
            
            NumRob = 500;
% %             rng(15)
% %             tempReg = randi(length(this.RegionID),1);
% %             temp1 = find(this.freespace(:,3) ==tempReg);
% %             temp2 = randi(numel(temp1),NumRob,1);
% %             this.loc = zeros(NumRob,2);
% %             this.loc(:,1:2) = this.freespace(temp1(temp2),1:2);
%             this.loc = zeros(NumRob,2);
%             this.loc(:,1:2) = this.freespace(randi((length(this.freespace)),NumRob,1),1:2);
%             figure
%             imshow(this.BW);
%             hold on
%             hRob = scatter(this.loc(:,2),this.loc(:,1),3,'filled');
%             hRob.CData = [1 0 0];
%             while 1
%                 Idx1 = knnsearch(this.Path(:,1:2),this.loc);
%                 [~, Idx2] = max(this.Path(Idx1,4));
%                 mov_tar = Idx2;
%                 p1 = this.loc(mov_tar,1:2);
%                 p2 = this.Path(Idx1(Idx2),1:2);
%                 p2_pat = this.Path(Idx1(Idx2),3);
% 
%                 while (p1==p2) ~=[1, 1]
%                     prev = this.loc;
%                     theta = atan2(p2(1)-p1(1),p2(2)-p1(2));
%                     this.loc(:,1:2) = this.loc(:,1:2) + [round(sin(theta)).*ones(NumRob,1),...
%                         round(cos(theta)).*ones(NumRob,1)];
%                     Idx =find(diag(this.BW(this.loc(:,1),this.loc(:,2)))==0);
%                     this.loc(Idx,:) = prev(Idx,:);
%                     p1 = this.loc(mov_tar,1:2);
%                     hRob.XData = this.loc(:,2);
%                     hRob.YData = this.loc(:,1);
%                     pause(0.01);
%                 end
% 
%                 while this.Path(p2_pat,3)>0
%                     p2 = this.Path(p2_pat,1:2);
%                     p2_pat = this.Path(p2_pat,3);
%                     prev = this.loc;
%                     ds =p2-p1;
%                     this.loc(:,1:2) = this.loc(:,1:2) + [ds(1).*ones(NumRob,1),...
%                         ds(2).*ones(NumRob,1)];
%                     Idx =find(diag(this.BW(this.loc(:,1),this.loc(:,2)))==0);
%                     this.loc(Idx,:) = prev(Idx,:);
%                     p1 = this.loc(mov_tar,1:2);
%                     hRob.XData = this.loc(:,2);
%                     hRob.YData = this.loc(:,1);
%                     pause(0.00001);
% 
%                 end
%             end
       end
       
       function LocalControl(this)
            temp = knnsearch(this.Path(:,1:2),this.localPath(:,1:2));

            temp1 = ismember(this.freespace(:,4),temp);
            this.localfsp = this.freespace(temp1,:);
% % %             figure
% % %             RGB = double(cat(3, ~this.BW, ~this.BW, ~this.BW));
% % %             RGB(:,:,1) = RGB(:,:,1).*182./255+ double(this.BW);
% % %             RGB(:,:,2) = RGB(:,:,2).*228./255+ double(this.BW);
% % %             RGB(:,:,3) = RGB(:,:,3).*255./255+ double(this.BW);
% % %             imshow(RGB)
% % %             hold on
% % %             h = scatter(this.localfsp(:,2),this.localfsp(:,1),2,'filled')
% % %             h.CData = [0.85,0.7,1];
            NumRob = 500;

            this.loc = zeros(NumRob,2);
            this.loc(:,1:2) = this.localfsp(randi((length(this.localfsp)),NumRob,1),1:2);
            
            localBW = (this.BW & 0);
            for ii = 1:length(this.localfsp)
                localBW(this.localfsp(ii,1),this.localfsp(ii,2)) = 1;
            end

            RGB = double(cat(3, ~localBW, ~localBW, ~localBW));
            RGB(:,:,1) = RGB(:,:,1).*182./255+ double(localBW);
            RGB(:,:,2) = RGB(:,:,2).*228./255+ double(localBW);
            RGB(:,:,3) = RGB(:,:,3).*255./255+ double(localBW);
%             figure
%             imshow(RGB);
%             hold on
            temp = find(this.localPath(:,4)<=100+1/3*this.channel_width);
            Idx = knnsearch(this.Path(:,1:2),this.localPath(temp,1:2));
            temp = ismember(this.freespace(:,4),Idx);
            
%             hshadow = scatter(this.freespace(temp,2),this.freespace(temp,1),10,'filled');
%             hshadow.CData = [0.85,0.7,1];
%             htar =scatter(this.tar(1,2),this.tar(:,1),100,'p','filled');
%             htar.CData = [1 0 0];
%             hRob = scatter(this.loc(:,2),this.loc(:,1),3,'filled');
%             hRob.CData = [0 0 0];
            totalcost = 1e10;
            Idx1 = knnsearch(this.localPath(:,1:2),this.loc);
            findmedial = 0;
            nstep = 0;
            
            tic
            while totalcost > NumRob*(100+1/3*this.channel_width)
                [~, Idx2] = max(this.localPath(Idx1,4));
                mov_obj = Idx2;
                p1 = this.loc(mov_obj,1:2);
                p2 = this.localPath(Idx1(Idx2),1:2);
                p2_pat = this.localPath(Idx1(Idx2),3);

                while findmedial ==0
                    prev = this.loc;
                    theta = atan2(p2(1)-p1(1),p2(2)-p1(2));
                    this.loc(:,1:2) = this.loc(:,1:2) + [round(sin(theta)).*ones(NumRob,1),...
                        round(cos(theta)).*ones(NumRob,1)];
                    Idx =find(diag(localBW(this.loc(:,1),this.loc(:,2)))==0);
                    this.loc(Idx,:) = prev(Idx,:);
                    p1 = this.loc(mov_obj,1:2);
%                     hRob.XData = this.loc(:,2);
%                     hRob.YData = this.loc(:,1);
                    nstep = nstep + 1;
                    pause(0.00001);
                    Idx1 = knnsearch(this.localPath(:,1:2),this.loc);
                    totalcost = sum(this.localPath(Idx1,4));
                    if (p1 == p2)==[1 1]
                        findmedial =1;
                    end
                end

                while this.localPath(p2_pat,3)>0
                    p2 = this.localPath(p2_pat,1:2);
                    p2_pat = this.localPath(p2_pat,3);
                    prev = this.loc;
                    ds =p2-p1;
                    this.loc(:,1:2) = this.loc(:,1:2) + [ds(1).*ones(NumRob,1),...
                        ds(2).*ones(NumRob,1)];
                    Idx =find(diag(localBW(this.loc(:,1),this.loc(:,2)))==0);
                    this.loc(Idx,:) = prev(Idx,:);
                    p1 = this.loc(mov_obj,1:2);
%                     hRob.XData = this.loc(:,2);
%                     hRob.YData = this.loc(:,1);
                    pause(0.0001);
                    Idx1 = knnsearch(this.localPath(:,1:2),this.loc);
                    totalcost = sum(this.localPath(Idx1,4));
                    nstep = nstep + 1;
                end
                Idx1 = knnsearch(this.localPath(:,1:2),this.loc);
                if ~isempty(this.basescore)
                    if nstep> 10*this.basescore 
                        break;
                    end
                end
            end
            
            toc
            disp(datetime('now'))
            if isempty(this.basescore)
                this.basescore = nstep;
            else
                this.Expt(1,this.col0,this.rank0) =this.tar(1);
                this.Expt(2,this.col0,this.rank0) =this.tar(2);
                this.Expt(this.row0,this.col0,this.rank0) = nstep;
                assignin('base','Expt',this.Expt)
                fprintf('Experiment %d, sample %d, region %d \n',this.row0-2,this.col0,this.rank0);

            end
            pause(1)
            close all
       end
       
   end


end