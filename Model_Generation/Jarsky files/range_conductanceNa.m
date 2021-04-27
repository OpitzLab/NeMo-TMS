function vec_gNa = range_conductanceNa(nainfo,tree,options)
% nainfo: gbar, slope

if nargin < 1 || isempty(nainfo)
    warning ('You need to define the information for the Na conductance disribution')
end

if numel(nainfo.region) == 1
    isregion = find(strcmp(tree.rnames,nainfo.region));
    isregion = tree.R == isregion; 
elseif numel(nainfo.region) > 1
    for r = 1:numel(nainfo.region)
        if ~contains(nainfo.region{r},'tuft')
            isregion(r,1) = find(strcmp(tree.rnames,nainfo.region{r}));
            isregions(:,r) = tree.R == isregion(r,1);
        end
    end
    isregion = sum(isregions,2);
else
    warning('Define the regions where the Na+ conductance will be distributed as a function of distance from the soma')
end

path_len = Pvec_tree(tree);
vec_gNa = (1:size(tree.X))'; 
fold_increase = 1.5;
max_pathlen = max(path_len.*isregion);
istuft = tree.R == find(strcmp(tree.rnames,'tuft'));

if contains(options,'-sE')          % Strong excitability model: increasing Na conductance with distance
    slope = (nainfo.gbar*fold_increase - nainfo.gbar)/(nainfo.gbar*max_pathlen);
    for node = 1:size(path_len)
        vec_gNa(node) = nainfo.gbar*(1+path_len(node)*slope);   % vector with the increasing values of hh Na conductance with distance
    end
    vec_gNa(istuft) = nainfo.gbar*fold_increase;  % tuft area gets the maximum conductance value
    
elseif contains(options,'-wE')      % Weak excitability model: uniform Na conductance (naslope = 0)
    for node = 1:size(vec_gNa)
        vec_gNa(node) = nainfo.gbar;
    end
end
vec_gNa(~isregion & ~istuft) = NaN;   % for the areas that are not part of the areas that will be ranged, no value is set
