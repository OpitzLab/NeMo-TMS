function tree = strip_axon(tree)
%Strips the axon away from the rest of the cell morphology.

axon_nodes = [];
axon_R = 2;
for i = 1:numel(tree.rnames)
    if strcmp(tree.rnames{i}, 'axon') || strcmp(tree.rnames{i}, 'Axon')
        axon_R = i;
    end
end

for i = 1:length(tree.R)
    if tree.R(i) == axon_R
        axon_nodes = [axon_nodes, i];
    end
end

tree = delete_tree(tree, axon_nodes);

%Rectify region assignments
for i = 1:length(tree.R)
    if(tree.R(i) == 3)
        tree.R(i) = 4;
    elseif tree.R(i) == 2
        tree.R(i) = 3;
    end
end

tree.rnames = {'soma', 'axon', 'basal', 'apical'};

end