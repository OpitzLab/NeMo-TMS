function tree = strip_noregion(tree)
%Strips the noregion thingy

nor_nodes = [];
nor_R = 5;


for i = 1:length(tree.R)
    if tree.R(i) == nor_R
        nor_nodes = [nor_nodes, i];
    end
end

tree = delete_tree(tree, nor_nodes);

%Rectify region assignments
for i = 1:length(tree.R)
    if(tree.R(i) == 1)
        tree.R(i) = 5;
    elseif tree.R(i) == 3
        tree.R(i) = 1;
    elseif tree.R(i) == 4
        tree.R(i) = 3;
    end
    
    if tree.R(i) == 5;
        tree.R(i) = 4;
    end
end

tree.rnames = {'soma', 'axon', 'basal', 'apical'};

end