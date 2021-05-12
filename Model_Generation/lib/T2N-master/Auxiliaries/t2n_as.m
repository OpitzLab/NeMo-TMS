function neuron = t2n_as(x,neuron)
% This function adds all missing definition fields to each neuron structure 
% in a neuron cell array by telling T2N to use the same fields as those of 
% the xth neuron structure in this array. Alternatively it fills up missing
% definitions in a single t2n structure. Useful if many neuron simulations 
% that should be run in parallel use partly the same definitions (e.g. same 
% mechanisms or recordings). Apply this function after everything that is
% different from the first neuron instance is defined.
%
% INPUTS
% x                 index to the neuron instance to use missing definitions from
% neuron            t2n neuron cell array of structures or single structure (see documentation)
%
% OUTPUTS
% neuron            t2n neuron structure or cell array of structures with
%                   filled up fields
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

fields = {'tree','mech','pp','con','record','play','APCount','params','tvecs','custom'};
if nargin < 1 || isempty(x)
    x = 1;
end

if nargin < 2
    for f = 1:numel(fields)
        if strcmp(fields{f},'tree')
            neuron.(fields{f}) = sprintf('sim%d',x);
        else
            neuron.(fields{f}) = x;
        end
    end
else
    if iscell(neuron)
        fields = intersect(fields,fieldnames(neuron{x}));
        for n = setdiff(1:numel(neuron),x)
            for f = 1:numel(fields)
                if strcmp(fields{f},'tree')
                    if ~isfield(neuron{n},fields{f}) || ischar(neuron{n}.(fields{f}))  % only replace value if not defined or pointing to another simulation
                        neuron{n}.(fields{f}) = sprintf('sim%d',x);
                    end
                elseif ~isfield(neuron{n},fields{f}) || isnumeric(neuron{n}.(fields{f}))  % only replace value if not defined or pointing to another simulation
                    neuron{n}.(fields{f}) = x;
                end
            end
        end
    elseif isstruct(neuron)
        fields = setdiff(fields,fieldnames(neuron));
        for f = 1:numel(fields)
            if strcmp(fields{f},'tree')
                neuron.(fields{f}) = sprintf('sim%d',x);
            else
                neuron.(fields{f}) = x;
            end
        end
    end
end
end
