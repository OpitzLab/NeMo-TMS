function [nodes,layersAP,plane,tree] = CA1pyramidalcell_divide(tree,addsoma,plot)
% Function to divide CA1 cells into differen regions and layers
% [nodes,layersAP,perp_dist,tree] 
% INPUT:
% Tree              tree structure
% Options           '-s'   plot all the different partitions

% OUTPUT: all logical arrays that specify the nodes in different regions
% soma_nodes        tree nodes from the somatic region
% basal_nodes       tree nodes from the basal region
% apical_nodes      tree nodes from the whole apical region
% ******************
% trunk_nodes       tree nodes from the main apical dendrite (until it
%                   reaches the s. lacunosum-moleculare layer
% oblique_nodes     tree nodes from the subtrees that stem from the trunk
% tuft_nodes        non-oblique (and non-trunk) apical tree nodes
% peritrunk_nodes   from Poirazi et al., 2003
% ******************
% SR_nodes          tree nodes that are above the SO-SP limit and below the 
%                   SR-SLM limit 
% SLM_nodes         tree nodes that are above the SR-SLM limit
% nodes_100         tree apical nodes that have less that 100um path length 
%                   or perpendicular length) 
% nodes_300         tree apical nodes that have less that 300um path length 
%                   or perpendicular length)
% nodes_350         tree apical nodes that have less that 350um path length 
%                   or perpendicular length)
% nodes_500         tree apical nodes that have less that 500um path length 
%                   or perpendicular length)
% nodes_rest        rest of the tree apical nodes (more than 500um path  
%                   length or perpendicular length)
% **************************************************************************

if (nargin < 2)||isempty(addsoma)
    addsoma = ''; % {DEFAULT: no option}
end
if (nargin < 3)||isempty(plot)
    plot = ''; % {DEFAULT: no option}
end

%% Resample the trees and clean the tree
% tree = resample_tree(tree,1);
% tree = repair_tree(tree);
tree = tran_tree(tree,1);   % Center the root to the 0,0,0 coordinates
%% If a pyramidal soma wants to be added:
if contains (addsoma, 'yes')
    tree = soma_tree(tree,15,25,'-r');
    % average pyramidal soma diameter = 15 um (from Hippocampus book)
end
%% Storage of the nodes of each region:
soma_nodes = tree.R ==find(strcmp(tree.rnames,'soma'));
basal_nodes = tree.R ==find(strcmp(tree.rnames,'basal'));
apical_nodes = tree.R ==find(strcmp(tree.rnames,'apical')); % doesn't include soma nodes
%% Calculate the 3D orthogonal distance regression line
% minimizes the sum of the euclidean distances of all points to the line
% Find the coordinates of the nodes I want to use for the regression: the
% proximal apical nodes of the tree -> 2/3 out of the total apical nodes
eucldist = eucl_tree(tree,1);
eucldist_ap = sort(eucldist(apical_nodes));
nodes_reg = ismember(eucldist, eucldist_ap(1:round((numel(eucldist_ap)/4)*3)));
counter = 1;
coordinates = ones(sum(nodes_reg),3);
for n = 1:numel(nodes_reg)
    if nodes_reg(n) == 1
        coordinates(counter,:) = [tree.X(n), tree.Y(n), tree.Z(n)];
        counter = counter + 1;
    end
end
% Force the regression line to pass through the root node coordinates
N = size(coordinates,1);    % number of nodes to align
root_point = [tree.X(1), tree.Y(1), tree.Z(1)]; % root = node 1
dcoor = bsxfun(@minus,coordinates,root_point);  % residuals (In the case of root = 0,0,0 it is the same as the matrix coordinates
C = (dcoor'*dcoor)/(N-1);           % variance-covariance matrix of X
[R,D] = svd(C,0);             % singular value decomposition of C; C=R*D*R'
% The direction of the best fit line corresponds to R(:,1) -> this is what
% I need to build the plane
% D(1,1) is the variance of dcoor after projection on R(:,1)
% Parametric equation of best fit line: L(t) = root_point + t*R(:,1)', where t is a real number
% Coefficient of determineation; R^2 = (explained variance)/(total variance)
D = diag(D);
R2 = D(1)/sum(D);
plane = R(:,1);
%% Calculate the plane perpendicular to the regression line 
% Plane: built usig the direction of the line and the root point
% Equation1: R(1,1)*(x-root_point(1)) + R(2,1)*(y-root_point(2)) +
% R(3,1)*(z-root_point(3)) = 0;
[x_plane1, y_plane1] = meshgrid(-100:1:100);  
z_plane1 = (- R(1,1)*(x_plane1 - root_point(1)) -  R(2,1)*(y_plane1 - root_point(2)))/R(3,1) + root_point(3);  
%% Calculate the apical dendritic length that contributes to s. radiatum
seglens     = len_tree(tree);
aplen       = seglens(apical_nodes);
aptreenodes = find(apical_nodes);
% Calculation of the perpendicular distance of all nodes to the plane SO-SP
for n = 1:numel(apical_nodes)
    Q(n,:) = [tree.X(n), tree.Y(n), tree.Z(n)]; % coordinates point of interest
    PQ(n,:) = Q(n,:) - root_point;
    % Dot product between line from PQ to P1 and normal of the plane:
    perp_dist(n,1) = PQ(n,:)*(R(:,1));
end
% hack around to get the rigth direction of the normal vector of the plane
% (for some reason it calculates the inverted one):
max_distalnode = find(eucldist == max(eucldist.*apical_nodes));
if perp_dist(max_distalnode) < 0    % if negative
    perp_dist = perp_dist*(-1);
end
m           = [perp_dist(apical_nodes),aptreenodes,aplen];   % sort the nodes using the euclidean distance
m_sort      = sortrows(m,1);
% Relative percentages of dendritic length for each apical region: 
% 67% apical dendritic length s. radiatum, and 33% s. lacunosum-moleculare 
% (taken from bibliography)
radlimit        = sum(aplen)*0.67;   
radnodes        = ones(numel(m_sort(:,1)),1);
lmnodes         = ones(numel(m_sort(:,1)),1);
previous_sum    = 0;
for r = 1:numel(m_sort(:,1))
    lensum = previous_sum + m_sort(r,3);
    if lensum <= radlimit
        radnodes(r) = m_sort(r,2);  % nodes that are part of s. radiatum 
    else
        lmnodes(r) = m_sort(r,2); % nodes that are part of s. lacunosum-moleulare
    end
    previous_sum = lensum;
end
% find the furthest node within the radnodes:
rad_nodes = radnodes(radnodes>1);
treenodes = (1:numel(tree.X))';
last_radnode = find(ismember(perp_dist,max(perp_dist(ismember(treenodes,rad_nodes)))));
coor_lastradnode = [tree.X(last_radnode), tree.Y(last_radnode), tree.Z(last_radnode)];
%% Find the plane that divided s. radiatum and s. lacunosum-moleculare
% This plane is parallel to the fisrt plane: substitute the root point for
% the new last coordinate point:
% Equation1: R(1,1)*(x-coor_lastradnode(1)) + R(2,1)*(y-coor_lastradnode(2)) +
% R(3,1)*(z-coor_lastradnode(3)) = 0;
[x_plane2, y_plane2] = meshgrid(-100:1:100);  
z_plane2 = (- R(1,1)*(x_plane2 - coor_lastradnode(1)) -  R(2,1)*(y_plane2 - coor_lastradnode(2)))/R(3,1) + coor_lastradnode(3);
%% Find the tree nodes above the plane SR-SLM
% Calculate again the perpendicular distance for the second plane:
for n = 1:numel(apical_nodes)
    if apical_nodes(n) == 1
        Q(n,:) = [tree.X(n), tree.Y(n), tree.Z(n)]; % coordinates point of interest
        PQ(n,:) = Q(n,:) - coor_lastradnode;
        % Dot product between line from PQ to P1 and normal of the plane:
        perp_dist1(n,1) = PQ(n,:)*(R(:,1));
    else
        perp_dist1(n,1) = 0;
    end
end
% hack around to get the rigth direction of the normal vector of the plane
% (for some reason it calculates the inverted one):
if perp_dist1(max_distalnode) < 0    % if negative
    perp_dist1 = perp_dist1*(-1);
end
nodes_SLM = perp_dist1 > 0;
%% Find the tuft branches
% Find the terminal tips in the SLM
term_points = T_tree(tree);
term_points_SLM = find(nodes_SLM.*term_points);
ipar = ipar_tree(tree);
for n = 1:numel(term_points_SLM)
    paths_term_points_SLM(n,:) = ipar(term_points_SLM(n),:);    % path length points of all terminal
end
% Find the nodes from the trunk
M = paths_term_points_SLM(:);
M = M(M>0);
U = unique(M);
trunknodes = U(2<histc(M,unique(M)));    % nodes shared from at least half of all SLM terminal nodes
% round(numel(term_points_SLM)/2)
% Trunk only goes till the edge of the SLM:
trunknodes = trunknodes(~ismember(trunknodes,find(nodes_SLM)));
trunk_nodes = (ismember(treenodes,trunknodes) - soma_nodes)>0;
% Some terminal points that "invade" the SLM could actually be oblique ->
% only a small portion of the branch invades the SLM.
% Condition: half of the total length of the branch should cross the layer border
counter = 0;
for branch = 1:size(paths_term_points_SLM,1)
    thisbranch          = paths_term_points_SLM(branch,:)'; % this includes trunk nodes
    thisbranch          = thisbranch(thisbranch(~ismember(thisbranch,trunknodes))>0);
    totlen(branch,1)    = sum(thisbranch>0);
    SLMnodes{branch}    = thisbranch(ismember(thisbranch,find(nodes_SLM)));
    SLMlen(branch,1)    = sum(SLMnodes{branch}>0);
    if (totlen(branch,1)/3)*2 < SLMlen(branch,1)    % if 2/3s of the branch are in the SLM
        tuftnodes(counter+1:numel(thisbranch)+counter,1) = thisbranch;
        counter = numel(tuftnodes);
    end
end
% Also, all the braches stemming from branching points in the SLM, are tuft branches 
SML_branchpoints = find(B_tree(tree).*nodes_SLM);
% choose the one closest to the soma:
SML_branchpoint_closest = find(perp_dist == min(perp_dist(SML_branchpoints)));
[addnodes,~] = sub_tree (tree, SML_branchpoint_closest);
SLMbranchnodes = find(addnodes);
tuftnodes = [SLMbranchnodes;tuftnodes];
tuftnodes = unique(tuftnodes);
tuft_nodes = (ismember(treenodes,tuftnodes) - trunk_nodes)>0;
oblique_nodes = (apical_nodes - trunk_nodes - tuft_nodes)>0;

%% Calculation of peritrunk nodes
trunk_branchpoints = find(B_tree(tree).*trunk_nodes);
idpar = idpar_tree(tree); % Get parent indices
for b = 1:length(trunk_branchpoints)
    children{b} = find(idpar==trunk_branchpoints(b));
    if ~ismember(children{b}(1),find(trunk_nodes))
        nodes_branch{b} = seglens.*(sub_tree(tree,children{b}(1)));
        prev_branch_len = 0;
        for n = 1:numel(nodes_branch{b})
            branch_len = nodes_branch{b}(n) + prev_branch_len;
            if branch_len > 50
                children_end{b} = n;
                break
            elseif n == numel(nodes_branch{b})
                children_end{b} = max(find(nodes_branch{b}>0));
            else
                prev_branch_len = branch_len;
            end
        end
        peritrunk_branch{b} = children{b}(1):1:children_end{b};
        
    elseif ~ismember(children{b}(2),find(trunk_nodes))
        nodes_branch{b} = seglens.*(sub_tree(tree,children{b}(2)));
        prev_branch_len = 0;
        for n = 1:numel(nodes_branch{b})
            branch_len = nodes_branch{b}(n) + prev_branch_len;
            if branch_len > 50
                children_end{b} = n;
                break
            elseif n == numel(nodes_branch{b})
                children_end{b} = max(find(nodes_branch{b}>0));
            else
                prev_branch_len = branch_len;
            end
        end
        peritrunk_branch{b} = children{b}(2):1:children_end{b};
    end
end
counter = 0;
for b = 1:numel(peritrunk_branch)
    peritrunk_branchnodes(counter+1:counter+numel(peritrunk_branch{b}),1) = peritrunk_branch{b};
    counter = numel(peritrunk_branchnodes);
end
peritrunk_nodes = ismember(treenodes,peritrunk_branchnodes);

%% Nodes contributing to the different layers: see methods
% 100um path length(or perpendicular length) =~ 10% length contribution to the apical tree
% 300um path length(or perpendicular length) =~ 47% length contribution to the apical tree
% 350um path length(or perpendicular length) =~ 10% length contribution to the apical tree
% 500um path length(or perpendicular length) =~ 20% length contribution to the apical tree
% rest path length(or perpendicular length)  =~ 13% length contribution to the apical tree
limit100        = sum(aplen)*0.08;
limit300        = sum(aplen)*0.48; 
limit350        = sum(aplen)*0.67; 
limit500        = sum(aplen)*0.85;
% now it makes sense to sort the nodes using path distance, so there will
% not be weird branches with the end part within a group and the beginning 
% within the next gourp 
pathlen         = Pvec_tree(tree);
m1              = [pathlen(apical_nodes),aptreenodes,aplen];   % sort the nodes using the euclidean distance
m1_sort         = sortrows(m1,1);
nodes100        = ones(numel(m1_sort(:,1)),1);
nodes300        = ones(numel(m1_sort(:,1)),1);
nodes350        = ones(numel(m1_sort(:,1)),1);
nodes500        = ones(numel(m1_sort(:,1)),1);
nodesrest       = ones(numel(m1_sort(:,1)),1);
previous_sum    = 0;
for r = 1:numel(m1_sort(:,1))
    lensum = previous_sum + m1_sort(r,3);
    if lensum <= limit100
        nodes100(r) = m1_sort(r,2);  % nodes that are part of s. radiatum 
    elseif lensum <= limit300
        nodes300(r) = m1_sort(r,2); % nodes that are part of s. lacunosum-moleulare
    elseif lensum <= limit350
        nodes350(r) = m1_sort(r,2);
    elseif lensum <= limit500
        nodes500(r) = m1_sort(r,2);
    else
        nodesrest(r) = m1_sort(r,2);
    end
    previous_sum = lensum;
end
nodes_100  = ismember(treenodes,nodes100(nodes100>1));
nodes_300  = ismember(treenodes,nodes300(nodes300>1));
nodes_350  = ismember(treenodes,nodes350(nodes350>1));
nodes_500  = ismember(treenodes,nodes500(nodes500>1));
nodes_rest = ismember(treenodes,nodesrest(nodesrest>1));
%% OUTPUT
layersAP = struct;
layersAP.n100  = nodes_100;
layersAP.n300  = nodes_300;
layersAP.n350  = nodes_350;
layersAP.n500  = nodes_500;
layersAP.nrest = nodes_rest;
nodes = struct;
nodes.soma = soma_nodes;
nodes.basal = basal_nodes;
nodes.apical = apical_nodes;
nodes.trunk = trunk_nodes;
nodes.peritrunk = peritrunk_nodes;
nodes.oblique = oblique_nodes;
nodes.tuft = tuft_nodes;
%% Visualization:
if contains(plot, '-s')
    x       = dcoor*R(:,1);    % project residuals on R(:,1) 
    x_min   = min(x);
    x_max   = max(x);
    dx      = x_max-x_min;
    Xa      = (x_min-0.05*dx)*R(:,1)' + root_point;
    Xb      = (x_max+0.05*dx)*R(:,1)' + root_point;
    X_end   = [Xa;Xb];
    figure;
    hold all;
    plot_tree(tree);
    plot3(X_end(:,1),X_end(:,2),X_end(:,3),'-r','LineWidth',1) % best fit line 
    mesh(x_plane1,y_plane1,z_plane1)
    mesh(x_plane2,y_plane2,z_plane2)
end

end
    

       