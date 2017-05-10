%% Map 0  A Swarm of 20 particles at random nodes of the tree
   
    % Particles location generation: pick a node in tree
    %swarm_gen = randi(size(tree,1),SwarmSize,1);
    swarm_gen = zeros(SwarmSize,1);
    ii = 1;
    while ii <= SwarmSize
        swarm_gen(ii) = randi(size(tree,1),1);
        if tree(swarm_gen(ii),4)>40
            ii = ii+1;
        end
    end
    ptc = zeros(SwarmSize,2);
    ptc(:,1:2) = tree(swarm_gen,1:2); 
    
