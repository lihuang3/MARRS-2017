%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% collision
%%   check to see if a node is in collsion with obstacles
function collision_flag = collision(node, parent, world)

collision_flag = 0;

if ((node(1)>world.NEcorner(1))...
    || (node(1)<world.SWcorner(1))...
    || (node(2)>world.NEcorner(2))...
    || (node(2)<world.SWcorner(2)))
  collision_flag = 1;
else
    px = [node(1),parent(1)];
    py = [node(2),parent(2)];
    if (~isempty(polyxpoly(world.ce(:), world.cn(:), px, py)))
        collision_flag =1;
    end
    
end


