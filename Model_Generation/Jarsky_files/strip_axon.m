function tree = strip_axon(tree)
%Strips the axon away from the rest of the cell morphology.

axon_nodes = [];
for i = 1:length(tree.R)
    if tree.R(i) == 2
        axon_nodes = [axon_nodes, i];
    end
end

tree = delete_tree(tree, axon_nodes);

for i = 1:length(tree.R)
    if(tree.R(i) == 3)
        tree.R(i) = 4;
    elseif tree.R(i) == 2
        tree.R(i) = 3;
    end
end

tree.rnames(4) = tree.rnames(3);
tree.rnames(3) = tree.rnames(2);

end