function vec_gKa = range_conductanceKa(kainfo,tree)
% nainfo: gka, slope, proximal_limit,

% proximal limit: calculation comes from the Jarsky et al., 2005 model
% the maximum path length of Jarsky tree is 818.3232um, therefore, 100um
% (which is the proximal limit they set) corresponds to the 12% of that
% maximum path length; and 500um, te limit corresponds to the 61% of max
% path length

if nargin < 1 || isempty(kainfo)
    warning ('You need to define the information for the A-type K+ conductance disribution')
end
if numel(kainfo.region) == 1
    isregion = find(strcmp(tree.rnames,kainfo.region));
    isregion = tree.R == isregion; 
elseif numel(kainfo.region) > 1
    for r = 1:numel(kainfo.region)
        if ~contains(kainfo.region{r},'tuft')
            isregion(r,1) = find(strcmp(tree.rnames,kainfo.region{r}));
            isregions(:,r) = tree.R == isregion(r,1);
        end
    end
    isregion = sum(isregions,2); 
else
    warning('Define the regions where the A-type K+ conductance will be distributed as a function of distance from the soma')
end
path_len = Pvec_tree(tree);
vec_gKa = struct;
vec_gKa.proximal = (1:size(tree.X))'; 
vec_gKa.distal = (1:size(tree.X))';
fold_increase = 6;
max_pathlen = max(path_len.*isregion);

for node = 1:size(path_len)
    slope = (kainfo.gka*fold_increase - kainfo.gka)/(kainfo.gka*max_pathlen);
    vecgKa(node,1) = kainfo.gka*(1+path_len(node)*slope);
end

% Get proximal and distal vectors:
isregion_100 = tree.R == find(strcmp(tree.rnames,'proxAp'));  % region within first 100 um
isregion_300 = tree.R == find(strcmp(tree.rnames,'middleAp'));  % region within first 300 um
isregion_500 = tree.R == find(strcmp(tree.rnames,'distalAp'));  % region within first 500 um
isregion_tuft = tree.R == find(strcmp(tree.rnames,'tuft'));
vec_gKa.proximal = vecgKa.*isregion_100;
vec_gKa.distal = vecgKa.*(isregion - isregion_100);
vec_gKa.distal(isregion_tuft) = kainfo.gka*fold_increase;

vec_gKa.proximal(~isregion) = NaN;
vec_gKa.distal(~isregion) = NaN;
