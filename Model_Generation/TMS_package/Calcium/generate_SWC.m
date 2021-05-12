%% Import data
index = importdata(['..' filesep '..' filesep 'Results' filesep 'NEURON' filesep 'locs' filesep 'index.dat']);
locs = importdata(['..' filesep '..' filesep 'Results' filesep 'NEURON' filesep 'locs' filesep 'locs_all_seg.txt']);
diams = index.data(:,1);

locations = locs(:,1:3);
parent_locs = locs(:,4:6);
%% find parent indices
parents = nan(size(parent_locs,1),1);
for i = 1:length(parents)
    if(isnan(parent_locs(i,1)))
        parents(i) = -1;
    else
        parents(i) = find(ismember(locations,parent_locs(i,:), 'rows'), 1);
    end
end
%% assign regions
regions = nan(length(parents),1);
for i = 1:length(parents)
    if contains(index.rowheaders(i), 'soma')
        regions(i) = 1;
    elseif contains(index.rowheaders(i), 'basal')
        regions(i) = 3;
    elseif contains(index.rowheaders(i), 'axon')
        regions(i) = 2;
    elseif contains(index.rowheaders(i), 'hill')
        regions(i) = 2;
    elseif contains(index.rowheaders(i), 'iseg')
        regions(i) = 2;
    elseif contains(index.rowheaders(i), 'node')
        regions(i) = 2;
    elseif contains(index.rowheaders(i), 'myelin')
        regions(i) = 2;
    else %presume it's apical otherwise
        regions(i) = 4;
    end
end
%% Initial SWC file
id = (1:size(locations,1))';
swc = [id regions locations diams/2 parents];
%% Adjust Order
% Some applications require the parent id to refer only to previously
% defined nodes.
% Here we fix the order first by adjusting the id and the parent id and 
% then sorting based on the adjusted id.
id_adj = nan(length(id),1);
parent_adj = nan(length(id),1);
id_adj(1) = 1;
parent_adj(1) = -1;
counter = 2;
while (sum(isnan(id_adj)) ~= 0)
    for ii = 1:length(id)
        % Add new node if not already defined but its parent is defined
        if isnan(id_adj(ii)) && ~isnan(id_adj(parents(ii)))
            id_adj(ii) = counter;
            parent_adj(ii) = id_adj(parents(ii));
            counter = counter + 1;
        end
    end
end

% sort the adjusted SWC file based on the adjusted id 
swc_adj = nan(size(swc));
[~,index] = sort(id_adj);
for ii = 1:length(index)
    swc_adj = [id_adj(index) swc(index,2:6)  parent_adj(index)];
end
%% Three-point soma representation
% Convert to 3-point soma if soma is single-point
% http://neuromorpho.org/SomaFormat.html
if sum(swc_adj(:,2) == 1) == 1
    swc_3p = [swc_adj(1,:); nan(2,7); swc_adj(2:end,:)];
    swc_3p(4:end,1) = swc_3p(4:end,1) + 2;
    idx = find(swc_3p(:,7)>1);
    swc_3p(idx,7) = swc_3p(idx,7) + 2;
    swc_3p(2:3,1) = 2:3;
    swc_3p(2:3,[2 7]) = [1 1; 1 1];
    swc_3p(2:3,3:6) = repmat(swc_3p(1,3:6),2,1);
    swc_3p(2:3,4) = swc_3p(2:3,4) + [-swc_3p(1,6); swc_3p(1,6)];
else
    swc_3p = swc_adj;
end
%% convert um to m
swc_3p(:,3:6) = swc_3p(:,3:6)*1e-6;
%% Export SWC file
output_file = ['..' filesep '..' filesep 'Results' filesep ...
    'Calcium' filesep 'neuron_out.swc'];
writematrix(swc_3p,output_file,'FileType','text','Delimiter','space');