%% Particle Swarm Generalization
%===============================

clear all
close all
clc

%%  Load RRT Map Data
load('data1226VeinsInBrain.mat')

%load('data1207.mat');   % Standard map
%load('dataT_owRRT1212.mat')     % T map
%load('data1129_2.mat')

% Including map, tree node, region division

rng(18)

% Swarm Population
Nindex = 7;
SwarmSize = 2^Nindex;

% Experiment Counter
Expt = 1;
NEP = 1;
% Experiment Data
ExptData=zeros(50,10);
%load('ExptVB1227.mat')
%%

if(Nindex<=10)
    for Expt = 1:NEP
    tic
    %% Initialization    
    SwarmSize = 2^Nindex;
    flagNP = 0; % flagNP = 1 <-- already enter the junction
    flagIjt = zeros(2,1);   
    % flagIjt(1) >=1  <-- switch to the global planner
    % flagIjt(1) < 1  <-- continue for the local planner
    flagIjt(1) = 1;
    % flagIjt(2) <-- save the most recent local planner info

    %% Plot map
    %open('fig1.fig'); % standard fig
    %open('fig2.fig'); % T fig
    open('figVB.fig');  % Veins in brain 
    run('map3.m')  % standard map
    %run('map0.m')     % T map
    
        fm = getframe;
        [img_temp, Xum] = frame2im(fm);
        imRow = size(img_temp,1);
        imCol = size(img_temp,2);
        
        simuproc1 = zeros(imRow,imCol,6000);
        simuproc2 = zeros(imRow,imCol,6000);
        simuproc3 = zeros(imRow,imCol,6000);
        
        
    figure(1); hold on
    fig1 = figure(1);
   % fig1.Position = [20 450 440 330];
    %title('Particle World');

    % Obtain particle locations.
    h = scatter(ptc(:,1),ptc(:,2),15,'filled');
    h.CData = [0,0.3,0.3];

    % h1 tracks the current object in solid red
    h1 = scatter(h.XData(1),h.YData(1),15,'filled');
    h1.CData = [1 0 0];

    % h2 indicates the temporary goal in hollow blue
%     h2 = scatter(h.XData(1),h.YData(1),35);
%     h2.CData = [0 0 1];
    
    
    if ExptData(Expt,Nindex) ==0
      


        %% Recursive Algorithm
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Algorithm Initialization
            discDiv = 5; % discount division --> every discDiv distance for a discount
            gamma = 0.3; % Discount factor
            Rjt = 100; % Junction rewards
            r0 = 0; % Rewards for unit distance

            Ijt_cnt = 0;    % region indicator with Ijt_cnt =1 the farthest region
            n = 1;     % algorithm step counter
            RegIntRange = 2;    % a threshold to switch to farthest-first heuristics
            TempSubReg = [];    % most recent visited branch in a region
            TSR_flag = 0;       % TSR_flag = 1 <-- stick to the most visited branch

            StepSize = 1;       % Step size
            NextParticle = 1;   % NextParticle =1 <-- the current goal has achieved,
                                % look for the next object to control

          % Total cost calculation
            pbinX =  1+floor(h.XData(:)/bin.unitX); 
            pbinY =  1+floor(h.YData(:)/bin.unitY);
            pbin = accumarray([pbinX,pbinY],1,[sizeb,sizeb]);
            pbin = (pbin.*binAvg);
            pbinCost = sum(pbin(:));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Np = 1;     % The current object in h plot (the scatter plot)

        while(pbinCost>10*SwarmSize)    % while the sum of cost is larger than 10 times of the swarm size 
            pause(0.001);
        fm = getframe;
        [img_temp, xum] = frame2im(fm);
        
        simuproc1(:,:,n) =img_temp(:,:,1);
        simuproc2(:,:,n) =img_temp(:,:,2);
        simuproc3(:,:,n) =img_temp(:,:,3);
            % Save the last step location
            Xpre = h.XData;
            Ypre = h.YData;

            if (NextParticle >= 1)  % Select new object to move
                Np = 0; % clear the current object

                %% Global planer

                % swarm distribution properties 
                SwarmProp = zeros(SwarmSize,4);
                % Col1: node number
                % Col2: Junction Region
                % Col3: branch number
                % Col4: expectated discounted rewards to the neareast junction

                Pdist = zeros(SwarmSize,1); % particle distribution -> in terms of region number

                % Particle distribution 
                for ii = 1:SwarmSize
                    new_p.coor(1:2) = [h.XData(ii),h.YData(ii)];    % the coordinate of obj ii 
                    new_p.range = 2.5;  
                    new_p.neighbor = 0;
                    new_p.binX =  1+floor(new_p.coor(1)/bin.unitX);
                    new_p.binY =  1+floor(new_p.coor(2)/bin.unitY);
                    search = 1;
                    % Locate the nearest tree node in RRT
                    run('NeighborSearch.m');
                    Pdist(ii) = tree(new_p.neighbor,8); 
                    SwarmProp(ii,1) = new_p.neighbor;
                    SwarmProp(ii,2) = Pdist(ii);
                end


                % Global planner and local planner
                while(~ (Ijt_cnt>0 && numel(find(Pdist == Ijt(Ijt_cnt))) > 0))
                    % WHILE the previous region is all clear
                    % IF flagIjt(1)>=1, switch to the global planner, and we look
                    % for the next region which still has robots
                    if flagIjt(1)>=1 || Ijt_cnt >=numel(Ijt)-RegIntRange
                        flagIjt(1) = 0; 
                        Ijt_cnt = 0;    % region number in Ijt <-- region array

                    % Start from the farthest region (Ijt[1])and if the region contains
                    % any robot then do aggregation
                        while (Ijt_cnt<numel(Ijt))
                            Ijt_cnt = Ijt_cnt+1;
                            if numel(find(Pdist == Ijt(Ijt_cnt))) > 0
                                break;
                            end

                        end
                    else
                    % IF flagIjt(i) < 1, stick to local planner --> clear the
                    % current region flagIjt(2)

                        Ijt_cnt = find(Ijt == flagIjt(2));
                        TSR_flag = 1;
                        TempSubReg = tree(CTID,9); % The last branch we were in
                        flagIjt(1) =  flagIjt(1)+1/3;   % incremental is a fraction
                        % rather than '1' because we want to keep tracking this
                        % team and send them to more inner region
                    end
                end



                % Particle clustering

                % SwarmProp(:,1) and SwarmProp(:,2) available
                % i.e. node number and region number available

                LocalSwm = [];

                % Info for the curent branch
                JTbranch_temp =zeros(size(JTbranch{Ijt_cnt},1),4);
                % COl1: number of members
                % Col2: dominent percentage
                % Col3: branch rewards

                % LocalSwm(:,1) saves the sequence # in SwarmProp, NOT nodes in the tree
                % SwarmProp(:,1) saves the node in the tree
                LocalSwm(:,1) = find(SwarmProp(:,2)==Ijt(Ijt_cnt));

                %% Local Planner
                for ii = 1:size(LocalSwm,1)

                    if tree(SwarmProp(LocalSwm(ii,1),1),6) ~=1   % find Major node
                        parent = tree(SwarmProp(LocalSwm(ii,1),1),5);
                    else
                        parent = SwarmProp(LocalSwm(ii,1),1);
                    end
                    parent0 = parent;

                   % branch member +1
                    JTbranch_temp( tree(parent0,9),1) = JTbranch_temp( tree(parent0,9),1)+1;
                    SwarmProp(LocalSwm(ii),3) = tree(parent0,9);
                    % Prior particle to move (inside junction region)
                    if tree(parent,7) >0 && (Ijt_cnt<numel(Ijt)-RegIntRange) %&& tree(parent,7)~=Ijt(end)     

                        Np = LocalSwm(ii);  % Identify the next object in the scatter plot
                        CTID = SwarmProp(LocalSwm(ii),1);   % nearest node to dock
                        break;
                    end




                    while tree(parent,7)==0 
                     % The parent of current node is not a junction node
                     % Trace back till find the nearest junction node.
                        parent = tree(parent,5); 
                    end

                    % Calculate expected reward of the particle 
                    SwarmProp(LocalSwm(ii),4) = r0*(tree(parent0,4)-tree(parent,4))...
                        +100*gamma.^((tree(parent0,4)-tree(parent,4))/discDiv);

                    % Calculate branch weight
                    JTbranch_temp( tree(parent0,9),3) = JTbranch_temp( tree(parent0,9),3)+ SwarmProp(LocalSwm(ii),4);


                end

                % Branch population percentage
                JTbranch_temp(:,2) = JTbranch_temp(:,1)./SwarmSize;

                if Np == 0 % No particles in the priority zone
                    if(TSR_flag ==1)
                        % Stick to the last visited branch
                        ActID1 = TempSubReg;
                        TSR_flag = 0;
                    else
                        % recalculate branch
                        [ActVal1 ActID1] = max(JTbranch_temp(:,3));
                    end
                    ActID2 = find(SwarmProp(:,3) ==ActID1);

                    % Incase the region is empty.
                    ActID1 = [];
                    while isempty(ActID2)
                       if isempty(ActID1)    
                            [ActVal1 ActID1] = max(JTbranch_temp(:,3)); 
                       else

                           JTbranch_temp(ActID1,3) =0;
                           [ActVal1 ActID1] = max(JTbranch_temp(:,3));
                       end
                       ActID2 = find(SwarmProp(:,3) ==ActID1);

                    end

                    % switch local planner (close region use farthest first)
                    if  (Ijt_cnt<numel(Ijt)-RegIntRange)
                        [ActVal3 ActID3] = max(SwarmProp(ActID2,4));
                        Np = ActID2(ActID3);    % # in scatter plot
                        CTID = SwarmProp(Np,1); % nearest node to dock   
                        parent0 = tree(CTID,5);
                    else

                        farthestCost = 0;
                        for ii = 1:numel(LocalSwm)

                            if tree(SwarmProp(LocalSwm(ii),1),4) > farthestCost
                                farthestCost = tree(SwarmProp(LocalSwm(ii),1),4);
                                Np = LocalSwm(ii);
                                CTID = SwarmProp(LocalSwm(ii),1);

                            end
                        end
                        farthestCost    
                        parent0 = tree(CTID,5);
                        if parent == 0
                           NextParticle =1;
                           Np = 0;
                           continue
                        end


                    end


                end


                theta = atan2((tree(CTID,2)-h.YData(Np)),(tree(CTID,1)-h.XData(Np)));
                StepSize = min(1,norm(tree(CTID,1:2)-[h.XData(Np),h.YData(Np)]));
                h.XData(:) = h.XData(:) + StepSize*cos(theta);
                h.YData(:) = h.YData(:) + StepSize*sin(theta);
                NextParticle = 0; 
                h.XData(Np) = tree(CTID,1);
                h.YData(Np) = tree(CTID,2);
                h1.XData =  h.XData(Np);
                h1.YData =  h.YData(Np);

                for ii= 1:SwarmSize
                    if ii ~=Np
                        if collision([h.XData(ii),h.YData(ii)],[Xpre(ii),Ypre(ii)],world) ==1
                           h.XData(ii) = Xpre(ii);
                           h.YData(ii) = Ypre(ii);
                        end
                    end
                end

                StepSize = 1;
            end


            Xpre = h.XData(:);
            Ypre = h.YData(:);


            n = n+1;
            if(mod(n,100)==0)
                disp('Total Step = ')
                disp(n)
            end


             if ( norm([h.XData(Np),h.YData(Np)]-tree(parent0,1:2)) < 1e-5)

                 if Ijt_cnt >= numel(Ijt)-RegIntRange
                    if(parent0>1)
                       CTID = parent0;
                       parent0 = tree(CTID,5);
                    else
                       NextParticle = 1; 
                    end
                 else
                    if tree(CTID,8) ~= tree(parent0,8)
                        flagIjt(2) = tree(parent0,8);
                        flagNP =1;
                        CTID = parent0;
                        parent0 = tree(CTID,5);
                    else
                        if flagNP == 1 
                            NextParticle = NextParticle + 1/2;
                            CTID = parent0;
                            parent0 = tree(CTID,5);

                            if NextParticle >= 1
                                flagNP =0;
                            end
                        else
                            CTID = parent0;
                            parent0 = tree(CTID,5);
                        end
                    end
                 end
             end

             if (CTID <= 1)
                disp('Particle Reach the Goal!');
                NextParticle = 1;
                continue

             else
                theta = atan2((tree(parent0,2)-h.YData(Np)),(tree(parent0,1)-h.XData(Np)));
                StepSize = min(1,norm(tree(parent0,1:2)-[h.XData(Np),h.YData(Np)]));
                h.XData(:) = h.XData(:) + StepSize*cos(theta);
                h.YData(:) = h.YData(:) + StepSize*sin(theta);

%                 h2.XData =  tree(parent0,1);
%                 h2.YData =  tree(parent0,2);
                h1.XData =  h.XData(Np);
                h1.YData =  h.YData(Np);
                for ii= 1:SwarmSize

                    if (ii ~=Np && collision([h.XData(ii),h.YData(ii)],[Xpre(ii),Ypre(ii)],world) ==1)
                         h.XData(ii) = Xpre(ii);
                         h.YData(ii) = Ypre(ii);
                    end
                end


                pbinX =  1+floor(h.XData(:)/bin.unitX); 
                pbinY =  1+floor(h.YData(:)/bin.unitY);


                pbin = accumarray([pbinX,pbinY],1,[sizeb,sizeb]);
                pbin = (pbin.*binAvg);

                pbinCost = sum(pbin(:));


             end


        end
        % 
         disp('Expt =')
         disp(Expt)
         disp('Complete! Total steps = ')
         disp(n)

         ExptData(Expt,Nindex) = n;

         disp('Swarm Size = ')
         disp(SwarmSize)


        pause(1.5);


        
    end
    figure(1);hold off
    close(figure(1))
    toc
    
    end
Nindex = Nindex+1;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of User Control Input