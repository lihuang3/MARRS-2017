%% Map 3  A Swarm of 16 particles at one region of the map
   
    % Particles location generation: pick a node in tree
    %RegionNode = find(tree(:,8)==Ijt(1));
    %swarm_gen = randi(numel(RegionNode),SwarmSize,1);
    swarm_gen = randi(TreeNode,SwarmSize,1);
    ptc = zeros(SwarmSize,2);
    %ptc(:,1:2) = tree(RegionNode(swarm_gen),1:2); 
    ptc(:,1:2) = tree((swarm_gen),1:2); 
