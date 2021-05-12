function t2n_showNetState(neuron,trees,out,tim)

for t = 1:numel(trees)
    if ~isfield(trees{t},'artificial') || ~trees{t}.artificial
        trees{t} = rot_tree(tran_tree(trees{t}),[],'-m3dY');
    else
        trees{t}.X = 0;
        trees{t}.Y = 0;
        trees{t}.Z = 0;
        trees{t}.dA = sparse(1);
    end
end
dd = spread_tree(trees);

col = colorme(numel(trees));

for t = 1:numel(trees)
    vec = t2n_get(out,'v',tim,'cell');
    
    if ~isfield(trees{t},'artificial') || ~trees{t}.artificial
        plot_tree(trees{t},vec,dd{t})
    else
        
    end
    
    
    % go through all defined point processes (pps)
    syn = fieldnames(neuron.pp{t});
    for s = 1:numel(syn)
        for n = 1:numel(neuron.pp{t}.(syn))
            % check which defined connections are handling this pp
            conCands = arrayfun(@(x) isfield(x.target,'cell') && x.target.cell==t && strcmp(x.target.pp,syn) && (~isfield(x.target,'ppg') || x.target.ppg==n),neuron.con);
            % get the cell source id for each connected pp
            synSources = cell2mat(arrayfun(@(x) ismember(neuron.pp{t}.(syn)(n).node(:),x.target.node)*x.source.cell,neuron.con(conCands),'uni',0));
            if any(sum(synSources~=0,2)>1)
                warning('At least one synapse seems to receive multiple connections. Only one of these connections is shown')
            end
            % get only one source cell index per pp (the maximum index)
            synSources = max(synSources{n},[],2);
            uSynSources = nonzeros(unique(synSources)); % ignore zeros (not connected when plotting synapses with the color of their source cell)
            for u = 1:numel(uSynSources)
                pointer_tree(trees{t},neuron.pp{t}.(syn)(n).nodes(synSources==uSynSources(u)),500,col{uSynSources(u)},dd{t},'-o')
            end
        end
    end
    
    
end



