a = find(tree(:,7)>0);
      b = find(abs(tree(a,1)-20.46)<0.1 & abs(tree(a,2)-99.83)<0.1)
      
      
      a(b)  = [];
      
h = scatter( tree(a,1),tree(a,2),50,'filled'); 

  targetH = scatter(99,5,50,'filled');
  targetH.CData = [1 0 0];