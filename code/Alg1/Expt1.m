%% Particle Swarm Generalization
%===============================

%#1 particles at tree nodes\
clear all
close all
clc

load('data1221Veins.mat')
%load('data1207.mat');   % Standard map
%load('dataT_owRRT1212.mat')     % T map
%load('data1129_2.mat')

% Including map, tree node, region division


% Swarm Population
rng(10)
Nindex = 6;
SwarmSize = 2^Nindex;

Expt = 1;
ExptData=zeros(30,10);
if(Nindex<=6)
    
SwarmSize = 2^Nindex;
flagNP = 0;

% Plot map
open('figVein.fig')
%open('fig1.fig'); % standard fig
%open('fig2.fig'); % T fig
run('map3.m')  % standard map
%run('map0.m')     % T map

figure(1); hold on
fig1 = figure(1);
fig1.Position = [20 450 440 330];
title('Particle World');

% Obtain particle locations.
h = scatter(ptc(:,1),ptc(:,2),20,'filled');
h.CData = [0,0.3,0.3];


h1 = scatter(h.XData(1),h.YData(1),20,'filled');
h1.CData = [1 0 0];
h2 = scatter(h.XData(1),h.YData(1),35);
h2.CData = [0 0 1];

%% End of Particle Swarm Generalization 
%

%% Algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
    discDiv = 5; % discount division --> every discDiv distance for a discount
    gamma = 0.3; % Discount factor
    Rjt = 100; % Junction rewards
    r0 = 0; % Rewards for unit distance
    
    Ijt_cnt = 0;
    n = 0;
    
    
    StepSize = 1;
    NextParticle = 1;
    
  % Entropy Computation
    grid = 32;   %refinement of map
    minvals = [0 0];%min(points);
    maxvals = [100, 100]; %max(points);
    rangevals = maxvals - minvals;
    pbinX =  1+floor(h.XData(:)/bin.unitX); 
    pbinY =  1+floor(h.YData(:)/bin.unitY);


    pbin = accumarray([pbinX,pbinY],1,[sizeb,sizeb]);
    pbin = (pbin.*binAvg);

    pbinCost = sum(pbin(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Np = 1;
while(pbinCost>10*SwarmSize)%(new_p.neighbor>0)
    pause(0.001);

    % Save the last step location
    Xpre = h.XData;
    Ypre = h.YData;
    
    if (NextParticle >= 1)  % Select new object to move
        Np = 0; % clear the target particle
        
        SwarmProp = zeros(SwarmSize,4);
        % Col1: node number
        % Col2: Junction Region
        % Col3: branch number
        % Col4: expectated rewards to the neareast junction
        
        Pdist = zeros(SwarmSize,1); % particle distribution -> in terms of region number
        
        % Particle distribution 
        %  o(???*o)
        for ii = 1:SwarmSize
            new_p.coor(1:2) = [h.XData(ii),h.YData(ii)];
            new_p.range = 2.5;
            new_p.neighbor = 0;
            new_p.binX =  1+floor(new_p.coor(1)/bin.unitX);
            new_p.binY =  1+floor(new_p.coor(2)/bin.unitY);
            search = 1;
            % Locate the nearest node
            run('NeighborSearch.m');
            Pdist(ii) = tree(new_p.neighbor,8);
            SwarmProp(ii,1) = new_p.neighbor;
            SwarmProp(ii,2) = Pdist(ii);
        end
        %  []~(???)~*
        
        if(~ (Ijt_cnt>0 && numel(find(Pdist == Ijt(Ijt_cnt))) > 0))
          Ijt_cnt = 0;    % region number in Ijt <-- region array
        
        % Start from the farthest region (Ijt[1])and if the region contains
        % more than 5% of the population then do aggregation
            while (Ijt_cnt<=numel(Ijt))
                Ijt_cnt = Ijt_cnt+1;
                if numel(find(Pdist == Ijt(Ijt_cnt))) > 0
                    break;
                end
            
            end
        
        end

        
        
        % Particle clustering
        %  o(???*o)
        
        % SwarmProp(:,1) and SwarmProp(:,2) available
        % i.e. node number and region number available
        
        LocalSwm = [];
        
        JTbranch_temp =zeros(size(JTbranch{Ijt_cnt},1),4);
        % COl1: number of members
        % Col2: dominent percentage
        % Col3: branch rewards
        
        % LocalSwm(:,1) saves the sequence # in SwarmProp, NOT nodes in the tree
        % SwarmProp(:,1) saves the node in the tree
        LocalSwm(:,1) = find(SwarmProp(:,2)==Ijt(Ijt_cnt));
        for ii = 1:size(LocalSwm,1)
            % !!!!!!!!!!!!!!!!!!
            if tree(SwarmProp(LocalSwm(ii,1),1),6) ~=1   % find Major node
                parent = tree(SwarmProp(LocalSwm(ii,1),1),5);
            else
                parent = SwarmProp(LocalSwm(ii,1),1);
            end
            parent0 = parent;
            
           % branch member +1
            JTbranch_temp( tree(parent0,9),1) = JTbranch_temp( tree(parent0,9),1)+1;
            SwarmProp(LocalSwm(ii),3) = tree(parent0,9);
            % Priority particle to move (inside junction region)
            if tree(parent,7) >0 && tree(parent,7)~=Ijt(end)
                
                Np = LocalSwm(ii);  % # in scatter plot
                CTID = SwarmProp(LocalSwm(ii),1);   % nearest node to dock
                break;
            end
            
            
            
             
            while tree(parent,7)==0 
             % The parent of current node is not a junction node
             % Trace back till find the nearest junction node.
                parent = tree(parent,5); 
            end
            
            % Compute expected reward of the particle 
            SwarmProp(LocalSwm(ii),4) = r0*(tree(parent0,4)-tree(parent,4))...
                +100*gamma.^((tree(parent0,4)-tree(parent,4))/discDiv);
            
            
            JTbranch_temp( tree(parent0,9),3) = JTbranch_temp( tree(parent0,9),3)+ SwarmProp(LocalSwm(ii),4);
            
%                 
% %             if tree(new_p.neighbor,4) > farthestCost
% %                 farthestCost = tree(new_p.neighbor,4);
% %                 Np = ii;
% %                 CTID = new_p.neighbor;
% % 
% %             end
        end
       %  []~(???)~*
       
        JTbranch_temp(:,2) = JTbranch_temp(:,1)./SwarmSize;
        
        if Np == 0 % No particles in the priority zone
            [ActVal1 ActID1] = max(JTbranch_temp(:,3));
            ActID2 = find(SwarmProp(:,3) ==ActID1);
            if(Ijt_cnt<numel(Ijt))    
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

                parent0 = tree(CTID,5);
                if parent == 0
                   NextParticle =1;
                   Np = 0;
                   continue
                end

                
            end

            
        end

%         parent = tree(CTID,5);
%         if parent == 0
%            NextParticle =1;
%            Np = 0;
%            continue
%         end

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


    if (1)
        n = n+1

%          % new_p.neighbor return the node near the centroid
%          new_p.coor(1:2) = [h.XData(Np),h.YData(Np)];
            
         if ( norm([h.XData(Np),h.YData(Np)]-tree(parent0,1:2)) < 1e-5)
                
             if Ijt_cnt == numel(Ijt)
                if(parent0>1)
                   CTID = parent0;
                   parent0 = tree(CTID,5);
                else
                   NextParticle = 1; 
                end
             else
                if tree(CTID,8) ~= tree(parent0,8)
                    flagNP = 1; 
                    CTID = parent0;
                    parent0 = tree(CTID,5);
                else
                    if flagNP == 1 
                        NextParticle = NextParticle + 1/40;
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

            h2.XData =  tree(parent0,1);
            h2.YData =  tree(parent0,2);
            h1.XData =  h.XData(Np);
            h1.YData =  h.YData(Np);
            for ii= 1:SwarmSize

                if (ii ~=Np && collision([h.XData(ii),h.YData(ii)],[Xpre(ii),Ypre(ii)],world) ==1)
                     h.XData(ii) = Xpre(ii);
                     h.YData(ii) = Ypre(ii);
                end
            end
            
%             dist = tree(parent0,1:2)-[h.XData(Np),h.YData(Np)]; 
%             dist = [h.XData(Np),h.YData(Np)]+dist./norm(dist).*StepSize; 
%             h.XData(Np) = dist(1);
%             h.YData(Np) = dist(2);




    %             xidx = 1 + round((h.XData(:) - minvals(1)) ./ rangevals(1) * (grid-1));
    %             yidx = 1 + round((h.YData(:) - minvals(2)) ./ rangevals(2) * (grid-1));
    %             density = accumarray([yidx, xidx], 1, [grid,grid]);  %note y is rows, x is cols
    % 
    %             entro = entropy(uint8(density));

    % 
    %             figure(3)
    %             scatter(n,entro,12,'filled');

            pbinX =  1+floor(h.XData(:)/bin.unitX); 
            pbinY =  1+floor(h.YData(:)/bin.unitY);

            
            pbin = accumarray([pbinX,pbinY],1,[sizeb,sizeb]);
            pbin = (pbin.*binAvg);
            %disp('Current Cost = ')
            pbinCost = sum(pbin(:));
    %             
    % 
    %             figure(6)
    %             scatter(n,pbinCost,12,'filled');
    %             


         end
 
     end
end
% 
 disp('Expt =')
 disp(Expt)
 disp('Complete! Total steps = ')
 disp(n)
 disp('Swarm Size = ')
 disp(SwarmSize)
 %ExptData(Expt,2) = numel(LocalSwm);
 ExptData(Expt,Nindex) = n;
 %ExptData(Expt,3) = 100*numel(LocalSwm)/SwarmSize;
% n = 1;
% disp('Swarm Size = ')
% disp(SwarmSize)

if (mod(Expt,30)==0 )
    Nindex = Nindex+1;
    Expt = 1;
else
    Expt = Expt +1;
end
pause(1.5);


figure(1);hold off
close(figure(1))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of User Control Input