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
      EBW
      Skel
      goal
      max_step
      ConnMat
      Pathway
      track
      mission_increment
      RegionID
      RegionVal
      trackL
      rtrack
   end
   
   methods
       %% class constructor
       function this = TestObj(map)
           
           %% Load Map
           testmap = strcat(map.name,'.fig');
           open(testmap);
           this.fig = getframe;
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
           figure
           imshow(this.Skel);
           hold on
           f1 = scatter(this.BrchPts(:,2),this.BrchPts(:,1),20,'filled');
           f1.CData = [1 0 0];
           f2 = scatter(this.EdPts(:,2),this.EdPts(:,1),20,'filled');
           f2.CData = [0 0 1];
           
           
           %% Preprocessing part I
           
           this.MapProcess1;
           this.MapProcess2;
           this.MapProcess3;
           this.MapProcess4;
           this.MapProcess5;
           this.MapProcess6;
            
       end
       
       %% Assign cost (from the goal)to nodes
       function MapProcess1(this)   
            this.Pathway = uint16(this.Skel);                   % 'Pathway' is a cost-to-go value map
            cost = 1e2;                                         % base cost
            this.Pathway(this.goal(1),this.goal(2))=cost;       % goal cost
            frtr = this.goal;                                   % BFS frontier
            this.sdEdPts = [];
            
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

                    cost = cost +1;
                    % Add new frontier nodes to the var 'frtr'
                    for ij = 1:numel(row)       
                        frtr = [frtr; row(ij)+frtr(1,1), col(ij)+frtr(1,2)];     
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

                            else
                               % The branch pt is within the threshold distance, and it
                               % will be regarded as an end pt.    
                               temp = [row(1)+temp(1),col(1)+temp(2)];     % Update the intermediate pt
                               this.Brch0(temp(1),temp(2))=0;
                               this.Ed0(temp(1),temp(2))=1;
                               this.Ed0(temp1(1),temp1(2))=0;

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
            
            figure
            imshow(~(~this.BW+this.Skel)); hold on
            
            
            %plot(BrchPts(:,2),BrchPts(:,1),'r.');
            f1 = scatter(this.BrchPts0(:,2),this.BrchPts0(:,1),20,'filled');
            f1.CData = [1 0 0];
            f2 = scatter(this.EdPts0(:,2),this.EdPts0(:,1),20,'filled');
            f2.CData = [0 0 1];
            if ~isempty(this.sdEdPts)
                f3 = scatter(this.sdEdPts(:,2),this.sdEdPts(:,1),5,'filled');
                f3.CData = [0,0,0.8];
            end

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
           this.trackL = this.track;
           
           this.rtrack = this.track;
           
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
            figure(3);             
            hold on
            temp1 = this.ConnMat(1,2*1:2*1+1);
            % Initialization of target-motion Demo
            h = scatter(temp1(2),temp1(1),80,'filled');
            h.CData = [1 0 0];

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
                            h.XData = temp(2);
                            h.YData = temp(1);
                            pause(0.005);
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
               intit = jj;
               while ~flag1
                  temp = this.trackL(jj,2*ii-1);  
                  if abs(temp-this.rtrack(1,2*ii-1))<= thrhd && this.trackL(jj,2*ii)>0
                      jj = jj+1;
                  else
                      flag1 = 1;
                  end
               end
               this.rtrack(2,2*ii-1) = jj-intit+1;
               this.rtrack(2,2*ii) = sum(this.trackL(intit:jj,2*ii));

            end
       
       end
       
       
   end
   
end