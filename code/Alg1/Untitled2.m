a = find(tree(:,6)>0);
      b = find(abs(tree(a,1)-20.46)<0.1 & abs(tree(a,2)-99.83)<0.1)
      
      
      a(b)  = [];
      rng(687);
for ii = 1:max(tree(:,8))
    b = find(tree(a,8)==ii);
    h = scatter( tree(a(b),1),tree(a(b),2),50,'filled');
    h.CData = randi(255,1,3)./255;
end
  targetH = scatter(99,5,50,'filled');
  targetH.CData = [1 0 0];