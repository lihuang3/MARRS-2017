% Junction branch properties
%  o(???*o)
JTbranch = cell(10,1);
tree(:,9) = 0;
% Col1: average angle
% COl2: number of members
% Col3: dominent percentage
% Col4: branch rewards
for ii_JTC = 1: numel(Ijt)
    % currentJT contains all major nodes in the current region (Ijt[ijt_cnt])
    currentJT = MajorParentFilt(find( tree(MajorParentFilt,8)== Ijt(ii_JTC) ));

    for ii_JTC_temp = 1:numel(currentJT)
        parent = currentJT(ii_JTC_temp); 
        while tree(parent,7)==0 
        % If the parent of current node is not a junction node
        % Trace back till find the nearest junction node.
           parent = tree(parent,5); 
        end

        % branch angle  -20 ~ 340
        ActAng = 180./pi.*atan2(tree(parent,2)-tree(currentJT(ii_JTC_temp),2), ...
            tree(parent,1)-tree(currentJT(ii_JTC_temp),1));
        if ActAng <= -20
            ActAng = ActAng+360;
        end

        % Clustering
        JT_flag = 0;    % clustering flag
        if ~isempty(JTbranch{ii_JTC})
            for kk = 1:size(JTbranch{ii_JTC},1)
                if  abs(ActAng-JTbranch{ii_JTC}(kk,1))<=20
                    % branch avg angle
                    JTbranch{ii_JTC}(kk,1) = (JTbranch{ii_JTC}(kk,1).*JTbranch{ii_JTC}(kk,2)+ActAng)./(JTbranch{ii_JTC}(kk,2)+1); 
                    JTbranch{ii_JTC}(kk,2) = JTbranch{ii_JTC}(kk,2)+1;  % # of members
                    tree(currentJT(ii_JTC_temp),9)=kk;      % mark major node
                    JT_flag = 1;
                    break;
                end
            end
        end
        if JT_flag == 0 
           JTbranch{ii_JTC}(end+1,1) = ActAng;
           JTbranch{ii_JTC}(end,2) = 1;
           tree(currentJT(ii_JTC_temp),9)= size(JTbranch{ii_JTC},1);
        end

    end  
    JTbranch{ii_JTC}(:,3) = JTbranch{ii_JTC}(:,2)./sum(JTbranch{ii_JTC}(:,2));

end
% []~(???)~*