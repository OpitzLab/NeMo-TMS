function tree = addnewaxon(tree, D1, D2)
%Strips the axon away from the rest of the cell morphology.
%%
axon_nodes = [];
L = 30;
for i = 1:length(tree.R)
    if tree.R(i) == 2
        axon_nodes = [axon_nodes, i];
    end
end

idpar = idpar_tree(tree);
pvec = Pvec_tree(tree);
ax_root = 1;
dir = dir_tree(tree,'-n');
dendritic_dir = [1 0 0];%nanmean(dir)/sqrt(sum(nanmean(dir).^2));
soma_end = -1;

for i = 1:length(axon_nodes)
    if tree.R(idpar(axon_nodes(i))) ~= 2
        ax_root = axon_nodes(i);
    end
end

i = 1;
while soma_end == -1
    if tree.R(i) == 1 && tree.R(i+1) ~= 1
        soma_end = i;
    end
    i = i + 1;
end
soma_end = ceil(soma_end/2);
        

%grab diameter at axon root
% D1 = tree.D(ax_root);
% D2 = D1;
%  
% 
% %grab axon diameter at 60um
% i = 1;
% while pvec(axon_nodes(i)) <= 60
%      D2 = tree.D(axon_nodes(i));
%      i = i + 1;
% end

%Nuke old axon
tree = delete_tree(tree, axon_nodes);

%Fix region assignment
for i = 1:length(tree.R)
    if(tree.R(i) == 3)
        tree.R(i) = 4;
    elseif tree.R(i) == 2
        tree.R(i) = 3;
    end
end
%%

%add new axon
origin = [tree.X(soma_end), tree.Y(soma_end), tree.Z(soma_end)];
sec1 = origin - L*dendritic_dir;
sec2 = sec1 - L*dendritic_dir;

%%
[tree, idpar(end+1)] = insert_tree(tree, [1, 2, sec1, D1, soma_end], 'nix');
tree = insert_tree(tree, [1, 2, sec2, D2, idpar(end)], 'nix');

tree.rnames(4) = tree.rnames(3);
tree.rnames(3) = tree.rnames(2);





%%
end