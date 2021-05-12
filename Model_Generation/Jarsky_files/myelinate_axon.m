function tree = myelinate_axon(tree)
%Adds simple myelin to an existing axon, with nodes at every axon 
%bifurcation as well as between internodal myelin segments.
%Terminal branches below 20um are left unmyelinated; diameters are set to
%uniform 1um for myelin segments, 0.8um for unmyelinated.
%Note: Regions must be set to SWC standard, and axon root must be set to
%soma or dendrite region.
%Nicholas Hananeia, May 2020

%Presumes that regions already assigned with
%1: soma, 2: axon, 3: basal, 4: apical
L_hill = 5;
L_iseg = 5;
L_myelin = 100;
L_node = 1;
D_myelin = 1.5;
D_hill = 0.8;
D_node = 0.8;
D_iseg = 0.5;
length_threshold = 20; %Terminals shorter than this will not be myelinated

%Add axon subregions
tree.rnames = [tree.rnames,'hill'];
tree.rnames = [tree.rnames,'iseg'];
tree.rnames = [tree.rnames,'myelin'];
tree.rnames = [tree.rnames,'node'];

%Various data bits I'll need
idpar = idpar_tree(tree);
child = child_tree(tree);
axon_nodes = find(tree.R == 2);
len = len_tree(tree);
Pvec = Pvec_tree(tree);

%Find most distal axon node
[~, init_node] = max(Pvec(axon_nodes));
init_node = axon_nodes(init_node);
current_node = init_node;
path_to_root = [];


%% TRAVERSE AXON TREE
%Find the main axon path. This will be used to assign the initial segment
%and hillock.
while tree.R(current_node) == 2
    path_to_root = [path_to_root, current_node];
    current_node = idpar(current_node);
end

%Find all terminal paths. This will traverse the whole tree and allow
%assignment of myelin.
terminal_nodes = [];
paths = {};
terminal_lengths = [];
%First, we find all terminal branches
for i = 1:length(axon_nodes)
    %if it's a terminal node, add to list
    if child(axon_nodes(i))==0
        terminal_nodes = [terminal_nodes, axon_nodes(i)];
    end      
end
%Now traverse from each terminal node to the base of the axon.
for i = 1:length(terminal_nodes)
    current_node = terminal_nodes(i);
    index = 1;
    paths{i} = [];
    while (tree.R(current_node) == 2 || tree.R(current_node) > 4)
        paths{i} = [paths{i}, current_node];
        current_node = idpar(current_node);
        index = index + 1;
    end
end


%% MYELINATION AND NODES
%I nsert a node at every branch point and myelinate everything else.
for i = 1:length(paths)
    cum_child = 0;
    terminal_lengths(i) = 0; %This counts the length of the terminal branches
    terminal_flag = 0; %Flips to on as soon as we find a bifurcation
    for j = 1:length(paths{i})
        if cum_child ~= child(paths{i}(j))%In this case, we are at a branch
            tree.R(paths{i}(j)) = 8; %Assign to node
            cum_child = child(paths{i}(j)); %Reset child value
            terminal_flag = 1;
        else
            tree.R(paths{i}(j)) = 7; %Assign to myelin
        end
        if terminal_flag == 0
            terminal_lengths(i) = terminal_lengths(i) + len(paths{i}(j));
        end
        cum_child = cum_child + 1;
    end
end

%Internodal myelin segments spaced based on length
for i = 1:length(paths)
    paths{i} = flip(paths{i}); %We want to count from the root outwards
end
axon_length = 0;
for i = 1:length(paths)
    for j = 1:length(paths{i})
        if tree.R(paths{i}(j)) == 8 %If there's already a node there is a bifurcation
            axon_length = 0;
        elseif axon_length < L_myelin
            tree.R(paths{i}(j)) = 7;
        else 
            tree.R(paths{i}(j)) = 8;
            axon_length = 0;
        end
        axon_length = axon_length + len(paths{i}(j));
    end
end

%% UNMYELINATED TERMINAL SEGMENTS
short_terminals = [];
for i = 1:length(terminal_nodes)
    if terminal_lengths(i) < length_threshold
        short_terminals = [short_terminals, i];
    end
end

for i = 1:length(short_terminals)
    %Counting outwards in again
    paths{short_terminals(i)} = flip(paths{short_terminals(i)});
    j = 1;
    while tree.R(paths{short_terminals(i)}(j)) ~= 8
        tree.R(paths{short_terminals(i)}(j)) = 2;
        j = j + 1;
    end
    paths{short_terminals(i)} = flip(paths{short_terminals(i)});
end
        
%% AIS AND HILLOCK
%Find the longest path, this will be our main branch.
path_to_root = [];
for i = 1:length(paths)
    if length(paths{i}) > length(path_to_root)
        path_to_root = paths{i};
    end
end
%Now, assign AIS, hillock, and main axon regular nodes. 
%This is done last on purpose.
axon_length = 0;
for i = 1:length(path_to_root)
    if axon_length <= L_hill
        tree.R(path_to_root(i)) = 5;
    elseif axon_length < (L_hill + L_iseg)
        tree.R(path_to_root(i)) = 6;
    elseif tree.R(path_to_root(i)) ~= 8
        tree.R(path_to_root(i)) = 7;
    end
    axon_length = axon_length + len(path_to_root(i));
end

%% SET DIAMETER OF AXON SEGMENTS
for i = 1:length(tree.R)
    if tree.R(i) == 5
        tree.D(i) = D_hill;
    elseif tree.R(i) == 6
        tree.D(i) = D_iseg;
    elseif tree.R(i) == 7
        tree.D(i) = D_myelin;
    elseif tree.R(i) == 8 || tree.R(i) == 2
        tree.D(i) = D_node;
    end
end
