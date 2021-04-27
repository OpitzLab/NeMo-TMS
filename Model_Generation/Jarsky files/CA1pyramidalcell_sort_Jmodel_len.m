function tree    = CA1pyramidalcell_sort_Jmodel_len (tree, options)

if nargin < 2 || isempty (options)
    % DEFAULT
    options      = '-j';
end

% first sorts the tree into basal and apical. It also saves now the nodes from the axon, if there is one 
if any(strcmp(tree.rnames,'soma'))     % Only if there is already a region defined as 'soma' funtion will start
    tree_regions = tree.R;
    soma_nodes = tree_regions == find(strcmp(tree.rnames,'soma'));   % Logical array the size of the number of nodes of the tree, with ones in the positions where the some is defined
    axon_nodes = tree_regions == find(strcmp(tree.rnames,'axon'));
    index_axon_nodes = find(tree.R == 2);
    index_soma_nodes = find(tree.R == 1);
    if contains(options,'-axon')   % means the cells already have an axon
        hill_nodes = tree_regions == find(strcmp(tree.rnames,'hill'));
        iseg_nodes = tree_regions == find(strcmp(tree.rnames,'iseg'));
        node_nodes = tree_regions == find(strcmp(tree.rnames,'node'));
        myelin_nodes = tree_regions == find(strcmp(tree.rnames,'myelin'));
    end
    tree.rnames(1,:) = [];    % Delete all the regions names
    %% ************** division for Jarsky model model
    if contains(options,'-j')
        seglens = len_tree(tree);     % store the length of each segment containing the tree
        treenodes = (1:numel(tree.X))';   
        idpar = idpar_tree(tree,'-0');    % Index of the direct parent node to each node from the tree
        child_branches_len = child_tree(tree,seglens);    % Accumulated values of all children nodes to each node from tree
        % Values are integrated using seglens to get the length of each child branch
        soma_plus_attached_nodes = ismember(idpar,index_soma_nodes);   % Find which nodes have soma nodes as parents nodes
        index_soma_plus_attached_nodes = find(soma_plus_attached_nodes);    % Index number of the 'soma_attached_nodes'
        % Removing the nodes that are part of the soma in the soma_attached_nodes
        for row = 1:size(index_soma_plus_attached_nodes,1)
            if find (index_soma_plus_attached_nodes(row,1)==index_soma_nodes(:,1))
                continue
            elseif find (index_soma_plus_attached_nodes(row,1)==index_axon_nodes(:,1))
                continue
            else
                index_soma_attached_nodes(row) = index_soma_plus_attached_nodes(row);
            end
        end
        index_soma_attached_nodes = find(tree.R == 4); %Only look at apical dendrite here
        index_soma_attached_nodes = index_soma_attached_nodes(find(index_soma_attached_nodes~=0));     % Get rid of the positions repeated values
        idpar = idpar(find(idpar~=0));
        soma_attached_nodes_child_branches_len = child_branches_len(idpar);    % Lengths downstream of all branches leaving soma
        [~,index] = max(soma_attached_nodes_child_branches_len(:));
        apical_start_point = index_soma_attached_nodes(index); % From all the nodes attached to the soma, I choose the one who had the maximum langth, and define it as apical
        apical_nodes = sub_tree(tree,apical_start_point); % Sub_tree puts 1s in the indices of the cild nodes of 'apical_start_point'
        basal_nodes = logical(ones(size(tree.X))-apical_nodes-soma_nodes-axon_nodes);       
% ******* Dividing the apical region depending on length proportions      
        pathlen = Pvec_tree(tree);
        len = len_tree(tree);
        treenodes = (1:numel(tree.X))';
        m = [pathlen,len,treenodes].*apical_nodes;
        m_sort = sortrows(m,1);
        len_sort = m_sort(:,2);

        perc100 = 3.14/100;    % from the original Jarsky
        perc300 = 36.27/100;     % from the original Jarsky
        perc500 = 68.90/100;    % from the original Jarsky

        totaplen = sum(len.*apical_nodes);
        len100 = totaplen*perc100;
        len300 = totaplen*perc300;
        len500 = totaplen*perc500;

        previous_sum = 0;
        for n = 1:numel(treenodes)
            totsum = previous_sum + len_sort(n);
            if totsum <= len100
                nodes100(n) = m_sort(n,3);
            elseif totsum <= len300
                nodes300(n) = m_sort(n,3);
            elseif totsum <= len500
                nodes500(n) = m_sort(n,3);
            else
                nodestuft(n) = m_sort(n,3);
            end
            previous_sum = totsum;
        end
        nodes100 = nodes100(nodes100>0);
        nodes300 = nodes300(nodes300>0);
        nodes500 = nodes500(nodes500>0);
        nodestuft = nodestuft(nodestuft>0);

        close_nodes = ismember(treenodes,nodes100);
        middle_nodes = ismember(treenodes,nodes300);
        far_nodes = ismember(treenodes,nodes500);
        tuft_nodes = ismember(treenodes,nodestuft);

% ******* Defining the new regions
        tree.rnames{1} =  'soma';
        tree.rnames = [tree.rnames,'basal'];
        tree.rnames = [tree.rnames,'proxAp'];
        tree.rnames = [tree.rnames,'middleAp'];
        tree.rnames = [tree.rnames,'distalAp'];
        tree.rnames = [tree.rnames,'tuft'];
        tree.rnames = [tree.rnames, 'axon'];
        tree.R(soma_nodes)= find(strcmp(tree.rnames,'soma'));
        tree.R(basal_nodes)= find(strcmp(tree.rnames,'basal'));
        tree.R(close_nodes)= find(strcmp(tree.rnames,'proxAp'));
        tree.R(middle_nodes)= find(strcmp(tree.rnames,'middleAp'));
        tree.R(far_nodes)= find(strcmp(tree.rnames,'distalAp'));
        tree.R(tuft_nodes)= find(strcmp(tree.rnames,'tuft'));
        tree.R(axon_nodes) = find(strcmp(tree.rnames, 'axon'));

        if contains(options,'-axon')   % means the cells already have an axon
            tree.rnames = [tree.rnames,'hill'];
            tree.rnames = [tree.rnames,'iseg'];
            tree.rnames = [tree.rnames,'node'];
            tree.rnames = [tree.rnames,'myelin'];
            tree.R(hill_nodes)= find(strcmp(tree.rnames,'hill'));
            tree.R(iseg_nodes)= find(strcmp(tree.rnames,'iseg'));
            tree.R(node_nodes)= find(strcmp(tree.rnames,'node'));
            tree.R(myelin_nodes)= find(strcmp(tree.rnames,'myelin'));
        end
    end
end


