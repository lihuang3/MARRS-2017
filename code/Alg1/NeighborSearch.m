while(new_p.neighbor==0)
    switch (search)
        case 0
            
            NBid = (find( (particles.binX(1:node)==new_p.binX) & (particles.binY(1:node) == new_p.binY) > 0));
            for i_NS = 1: numel(NBid)
                 dist = norm( particles.coor( NBid(i_NS),1:2 )- new_p.coor );
                 if dist < new_p.range
                    new_p.neighbor = NBid(i_NS);
                    new_p.range = dist;
                 end
            end
            
        case 1
            binNextY =0;
            if ( (new_p.binY == 1+round(new_p.coor(2)/bin.unitY) ) && (new_p.binY-1 > 0))
               binNextY = -1;
            else
                if ((new_p.binY == round(new_p.coor(2)/bin.unitY) ) && (new_p.binY+1 <= sizeb))
                    binNextY = 1;
                end 
            end

            binNextX =0;
            if ( (new_p.binX == 1+round(new_p.coor(1)/bin.unitX) ) && (new_p.binX-1 > 0))
               binNextX = -1;
            else
                if ((new_p.binX == round(new_p.coor(1)/bin.unitX) ) && (new_p.binX+1 <= sizeb))
                    binNextX = 1;
                end 
            end

            for k_NS = 0:1  %% binY
                if(k_NS ==0 || binNextY ~= 0)
                    for j_NS = 0:1  %% binX
                        if(j_NS==0 || binNextX ~=0)
                            binX = new_p.binX + j_NS*binNextX;
                            binY = new_p.binY + k_NS*binNextY;
                            NBid = (find( (particles.binX(1:node)==binX) & (particles.binY(1:node) == binY) > 0));
                            for i_NS = 1: numel(NBid)
                                 dist = norm( particles.coor( NBid(i_NS),1:2 )- new_p.coor );
                                 if dist < new_p.range
                                    new_p.neighbor = NBid(i_NS);
                                    new_p.range = dist;
                                 end
                            end
                        end

                    end
                end
            end
            
            
        case 2
            for k_NS = -1:1 % search for neighbors in the bin and up and down bins
                if ( k_NS+new_p.binY > 0 && k_NS+new_p.binY <= sizeb )
                    for j_NS = -1:1 % search for neighbors in the bin and left/right bins

                        if ( j_NS+new_p.binX > 0 && j_NS+new_p.binX <= sizeb ) 
                            % search for all neighbors in the same bin
                            binX = new_p.binX+j_NS;
                            binY = new_p.binY+k_NS;
                            NBid = (find( (particles.binX(1:node)==binX) & (particles.binY(1:node) == binY) > 0));
                            for i_NS = 1: numel(NBid)
                                 dist = norm( particles.coor( NBid(i_NS),1:2 )- new_p.coor );
                                 if dist < new_p.range
                                    new_p.neighbor = NBid(i_NS);
                                    new_p.range = dist;
                                 end
                            end
                        end
                    end
                end
            end
           
        case 3
           for k_NS = -2:2 % search for neighbors in the bin and up and down bins
                if ( k_NS+new_p.binY > 0 && k_NS+new_p.binY <= sizeb )
                    for j_NS = -2:2 % search for neighbors in the bin and left/right bins

                        if ( j_NS+new_p.binX > 0 && j_NS+new_p.binX <= sizeb ) 
                            % search for all neighbors in the same bin
                            binX = new_p.binX+j_NS;
                            binY = new_p.binY+k_NS;
                            NBid = (find ((particles.binX(1:node)==binX) & (particles.binY(1:node) == binY) > 0));

                            for i_NS = 1: numel(NBid)
                                 dist = norm( particles.coor( NBid(i_NS),1:2 )- new_p.coor );
                                 if dist < new_p.range
                                    new_p.neighbor = NBid(i_NS);
                                    new_p.range = dist;
                                 end
                            end
                        end
                    end
                end
            end
            
        case 4

            for i_NS = 1: node
                dist =  norm(particles.coor( i_NS,1:2 )- new_p.coor) ;
                
                if dist < new_p.range
                      new_p.neighbor = i_NS;
                      new_p.range = dist;
                 end
                
            end
            
        
    end 
    search = search + 1;
end