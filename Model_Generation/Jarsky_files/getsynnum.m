function [nsyn, syn_ids, syndens] = getsynnum(tree,syn_dens,regions)
rng('shuffle') 
for r = 1:numel(regions)    
    ind = tree.R==find(strcmp(tree.rnames,regions{r}));
    nodes{r} = find(ind);
    len = len_tree(tree);
    region_len(r) = sum(len(ind));
    tot_nsyn = round(syn_dens*region_len(r));
    if tot_nsyn > numel(nodes{r})
        extra_nodes = tot_nsyn - numel(nodes{r});
        synids1 = datasample(nodes{r},(tot_nsyn-extra_nodes),'Replace',false);
        synids2 = datasample(nodes{r},extra_nodes,'Replace',false);
        synids{r} = [synids1;synids2];
    else
        synids{r} = datasample(nodes{r},tot_nsyn,'Replace',false);
    end
    syndens(r) = numel(synids{r})/region_len(r);
end
syndens = mean(syndens);
syn_ids = cell2mat(synids(:));
all_nodes = cell2mat(nodes(:));
not_synids = find(~ismember(all_nodes,syn_ids));
if syndens > (syn_dens + 0.0002)
    while syndens > (syn_dens + 0.0002)
        syn_ids = datasample(syn_ids,(numel(syn_ids)-1),'Replace',false);
        syndens = numel(syn_ids)/(sum(region_len));
    end
elseif syndens < (syn_dens - 0.0002)
    while syndens < (syn_dens - 0.0002)
        toadd = datasample(not_synids,1,'Replace',false);
        syn_ids = [syn_ids;toadd];
        syndens = numel(syn_ids)/(sum(region_len));
    end
end
nsyn = ismember((1:numel(tree.X))',syn_ids);
end
    
    
    