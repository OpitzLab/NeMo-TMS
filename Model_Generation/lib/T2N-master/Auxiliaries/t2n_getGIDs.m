function  [GIDs,neuron,mindelay] = t2n_getGIDs(neuron,tree,thesetrees)
% This function prepares the neuron structure to be used for the parallel
% NEURON environment. Also it returns the global IDs, which are necessary
% for the parallel environment.
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms (see documentation)
% tree              tree cell array with morphologies (see documentation)
% thesetrees        (optional) if not all trees of 'tree' are used, this is
%                   the index array to the used ones
%
% OUTPUTS
% GIDs              structure with information on all needed global IDs 
%                   that are initialized in parallel NEURON
% neuron            complemented neuron structure ready for parallel NEURON
% mindelay          minimum delay and time step of parallel network
%                   (calculated from minimal delay between synapses)
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if isfield(neuron.params,'dt') && (~isfield(neuron.params,'cvode') || ~neuron.params.cvode)
    mindelay = min(neuron.params.dt*2,10);
else
    mindelay = 10;  % default minimum delay and time step of parallel network
end
if ~exist('thesetrees','var') || isempty(thesetrees)
    thesetrees = 1:numel(tree);
end
% give each existing cell a unique GID and initialize some basic parameters
% this is necessary that gid_exist works
for t = 1:numel(thesetrees)
    if isfield(tree{thesetrees(t)},'artificial')
        GIDs(t) = struct('gid',t-1,'pp',[],'ppg',0,'watch','on','node',[],'cell',t,'threshold',[]);
    else
        GIDs(t) = struct('gid',t-1,'pp',[],'ppg',0,'watch','v','node',1,'cell',t,'threshold',[]);
    end
end
counter = 1+t; 

if isfield(neuron,'con')
    for c = 1:numel(neuron.con)
        if isfield(neuron.con(c).source,'pp') && ~isempty(neuron.con(c).source.pp)
            
            % add the point process which is used as a netcon as a new GID            
            ind = find(arrayfun(@(x) strcmp(x.pp,neuron.con(c).source.pp) & x.ppg == neuron.con(c).source.ppg & x.cell  == neuron.con(c).source.cell & isequal(x.node,neuron.con(c).source.node) & strcmp(x.watch,neuron.con(c).source.watch) ,GIDs),1,'first');
            if ~isempty(ind)
                neuron.con(c).source.gid = GIDs(ind).gid;
            else
                neuron.con(c).source.gid = counter-1;
                tmp = neuron.con(c).source;
                if isfield(neuron.con(c),'threshold')
                    tmp.threshold = neuron.con(c).threshold;
                else
                    tmp.threshold = [];
                end
                GIDs(counter) = tmp;
                counter = counter +1;
            end
        else
            % if not pp, check if GID of cell exist, check if it is
            % about the same node, check if the same variable is watched
            ind = arrayfun(@(x) isempty(x.pp) & x.cell  == neuron.con(c).source.cell & isequal(x.node,neuron.con(c).source.node) & strcmp(x.watch,neuron.con(c).source.watch) ,GIDs);
            if any(ind)
                neuron.con(c).source.gid = GIDs(ind).gid;
                GIDs(ind).threshold = neuron.con(c).threshold;
            else
                neuron.con(c).source.gid = counter-1;
                tmp = neuron.con(c).source;
                if isfield(neuron.con(c),'threshold')
                    tmp.threshold = neuron.con(c).threshold;
                else
                    tmp.threshold = [];
                end
                GIDs(counter) = tmp;
                counter = counter +1;
            end
        end
        % check if there is any netcon delay smaller than mindelay    
        if isfield(neuron.con(c),'delay')
            mindelay = min(mindelay,neuron.con(c).delay);
        end
    end
end
if mindelay <= 0
    error('The minimum delay of NetCons in parallel NEURON has to be greater than zero!')
end
[~,ind] = sort(cat(1,GIDs.cell));
GIDs = GIDs(ind);

end

