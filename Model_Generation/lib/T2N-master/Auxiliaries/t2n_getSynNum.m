function [nsyn,synids] = t2n_getSynNum(tree,syn_dens,regions)
% This function calculates the number and location of synapses at a given
% morphology given a synaptic density 'syn_dens'. The output synids can 
% then be used to place location-specific synapse mechanisms such as Exp2Syn.
%
% INPUTS
% tree              one morphology structure (see documentation)
% syn_dens          single scalar with the desired synaptic density 
%                   [#/length unit of morphologies]. If multiple regions
%                   are given (see next input), this can also be a vector
%                   with one entry for each region.
% regions           (optional) regions for which the synapse number should be
%                   calculated for. These can be anything from the 'rnames'
%                   field of the tree structure
%
% OUTPUTS
% nsyn              number of synapses at each node index of 'tree'
% synids            index to the node for each synapse to be implemented
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 3
    regions = tree.rnames;
end
if ~iscell(regions)
    regions = {regions};
end
if nargin < 2
    syn_dens = 1;
end
if numel(syn_dens) == 1 && numel(regions) > 1
    syn_dens = repmat(syn_dens,numel(regions),1);
end
if iscell(tree)
    if numel(tree)>1
        warning('Multiple trees as input detected. Using first one, ignoring others')
    end
    tree = tree{1};
end
nsyn = zeros(numel(tree.X),1);
synids = [];
for r = 1:numel(regions)
    
    ind = tree.R==find(strcmp(tree.rnames,regions{r}));
    
    sect = dissect_tree(tree);
    sect = sect(ind(sect(:,2)),:);
    ipar = ipar_tree(tree,[],sect(:,2));
    for n = 1:size(ipar,1)
        ipar(n,find(ipar(n,:)==sect(n,1)):end) = 1;  % make branch point and following nodes to the root refer to node 1 (root)
    end
    ipar = ipar(:,any(ipar ~= 1,1));  % delete unneeded columns
    len = len_tree(tree);        % get lengths of node segments
    synm = len(ipar)*syn_dens;   % calculate matrix with number of synapses per node
    synm(synm==0) = NaN;
    synm = round(synm*10)/10;    % round number of synapses
    flag = false;   % check flag if cumsum has been performed
    while 1 %any(synm(:)>=1)
        synids = cat(1,synids,ipar(synm>=1));
        nsyn(ipar(synm>=1)) = nsyn(ipar(synm>=1)) + 1;
        synm(synm>=1) = synm(synm>=1) -1;
        if ~any(synm(:)>=1)   % this has to be done to make it able to have non-integer spine densities
            if ~flag  
                for n = 2:size(synm,2)  % go through syn matrix and add up synapse numbers each time until value exceeds 1
                    add = synm(:,n-1);
                    add(add >= 1) = add(add >= 1) -1;
                    synm(:,n) = round((synm(:,n) + add)*10)/10;
                end
                flag = true;
            else
                break
            end
        end
    end
end