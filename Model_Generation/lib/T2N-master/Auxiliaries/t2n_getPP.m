function [varargout] = t2n_getPP(neuron,tree,pp)
% This function plots the locations of a given NEURON point process onto each
% morphology and optionally returns these values, too;
% INPUTS
% neuron    t2n neuron structure or cell array of structure (see
%           documentation)
% tree      TREES toolbox morphologies
% pp       string or cell array of strings with valid NEURON point process 
%           name, such as 'Exp2Syn'
%
% OUTPUT
% varVec    m-by-n cell array for m neuron specifications and n trees. each
%           contains a nxp logical array (n tree nodes, p = size(pp)) 
%           being true when the pp is defined at that node
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 3 || isempty(pp)
    pp = 'Exp2Syn';
    disp('Exp2Syn is used as point process')
end
if nargin < 1 || isempty(neuron)
    error 'neuron has not been defined or is no struct'
end
if nargin < 2 || isempty(tree)
    tree = {sample_tree};
    disp('Example tree is used')
end
if isstruct(tree)
    tree = {tree};
end
if isstruct(neuron)
    neuron = {neuron};
end
% add passive mechanism label to cm/Ra
if ~ iscell(pp)
    pp = {pp};
end
cols = colorme(numel(pp)+1);
varVec = cell(numel(neuron),numel(tree));

for n = 1:numel(neuron) % go through all neuron definitions
    for t = 1:numel(tree)  % go through all morphologies
        varVec{n,t} = false(numel(tree{t}.X),numel(pp));  % initialize the vector
        if isfield(neuron{n},'pp')   % check for mechanism definition
            ppnames = fieldnames(neuron{n}.pp{t}); % get ppnames
            for r = 1:numel(pp)  % go through point processes
                if any(strcmp(ppnames,pp))
                    varVec{n,t}(cat(1,neuron{n}.pp{t}.(ppnames{r}).node),r) = true;
                end
            end
        end
        if nargout == 0  % map locations on tree
            figure
            hold all
            if isfield(tree{t},'name')
                tname = tree{t}.name;
            else
                tname = sprintf('tree %d',t);
            end
            title(sprintf('Neuron Simulation: %d, tree: "%s", variable(s) "%s"',n,strrep(tname,'_','\_'),strrep(strjoin(pp,'|'),'_','\_')))
            plot_tree(tree{t});
            
            for p = 1:numel(pp)
                HP = pointer_tree(tree{t},varVec{n,t}(:,p), [], cols{p+1}, [], '-s');
                h(p) = HP(1);
            end
            legend(h,pp,'Location','best')
        end
    end
end
if nargout > 0
    varargout{1} = varVec;
end
end

