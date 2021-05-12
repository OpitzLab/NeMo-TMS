function  neuron = t2n_implicitSpines(neuron,scaleAll,tree)
% This function goes through all regions and scales the conductances (and 
% cm) according to the scaling information in region.spines. 
% Then it deletes the spines mechanism.
% 
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms
% scaleAll          Boolean determining if all conductances or only g_pas 
%                   and cm should be scaled
% tree              tree cell array with morphologies. Only necessary for
%                   ranged spines
%
% OUTPUTS
% neuron            modified t2n neuron structure
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 3
    tree = [];
end
if ~exist('scaleAll','var') || isempty(scaleAll)
    scaleAll = false;
end
if ~iscell(neuron)
    neuron = {neuron};
    cellflag = true;
else
    cellflag = false;
end

for n = 1:numel(neuron)
    for t = 1:numel(neuron{n}.mech)
        if ~isempty(neuron{n}.mech{t})
            fields = fieldnames(neuron{n}.mech{t});
            for f1 = 1:numel(fields)
                if isfield(neuron{n}.mech{t}.(fields{f1}),'spines')
                    scale = neuron{n}.mech{t}.(fields{f1}).spines.scale;
                    neuron{n}.mech{t}.(fields{f1}) = rmfield(neuron{n}.mech{t}.(fields{f1}),'spines');  %remove spine scale field from neuron structure
                    fields2 = fieldnames(neuron{n}.mech{t}.(fields{f1}));
                    if strcmp('range',fields{f1})  %thats complicated...go through all regions and find all gbars and make range variables out of it...
                        if isempty(tree)
                            error('Corresponding tree structure is needed as third input if spines should be scaled in a range instead of region manner!')
                        end
                        sscale = scale;
                        sscale(isnan(sscale)) = 1;  % in order to not make defined variables NaN
                        for ff2 = 1:numel(fields2) % go through variables defined as range variables
                            if strcmp(fields2{ff2},'pas')
                                neuron{n}.mech{t}.range.pas.g = neuron{n}.mech{t}.(fields{f1}).pas.g * sscale;
                                neuron{n}.mech{t}.range.pas.cm = neuron{n}.mech{t}.range.pas.cm * sscale;
                            elseif scaleAll
                                fields3 = fieldnames(neuron{n}.mech{t}.range.(fields2{ff2}));
                                for ff3 = 1:numel(fields3)
                                    if ~isempty(strfind(fields3{ff3},'bar'))
                                        neuron{n}.mech{t}.range.(fields2{ff2}).(fields3{ff3}) = neuron{n}.mech{t}.range.(fields2{ff2}).(fields3{ff3}) * sscale;
                                    end
                                end
                            end
                        end
                        rangefields = tree{t}.rnames(unique(tree{t}.R(~isnan(scale))));  % get all region fields defined in scale
                        for ff1 = 1:numel(rangefields)
                            ind = tree{t}.R == find(strcmp(tree{t}.rnames,rangefields{ff1}));   % get index to nodes that are in this region
                            fields2 = fieldnames(neuron{n}.mech{t}.(rangefields{ff1}));
                            
                            fields_all = fieldnames(neuron{n}.mech{t}.all);
                            fields_all = fields_all(cellfun(@(x) isempty(strfind(x,'_ion')),fields_all)); % sort out ion mechanisms
                            if ~isempty(setdiff(fields_all,fields2)) % if there are mechanisms defined in "all" regions that have not been addressed in the current region
                                if any(strcmp(fields_all,'spines'))
                                    warning('CAUTION! Field "spines" has been defined in region %s and in region "all"! Please check!',rangefields{ff1})
                                else
                                    morefields = setdiff(fields_all,fields2);
                                    for f2 = 1:numel(morefields)
                                        neuron{n}.mech{t}.(rangefields{ff1}).(morefields{f2}) = neuron{n}.mech{t}.all.(morefields{f2});  % property is now explicitly defined for that region so that it can be modified by spines information
                                    end
                                    fields2 = cat(1,fields2,morefields);
                                end
                                
                            end
                            
                            
                            for ff2 = 1:numel(fields2)
                                if strcmp(fields2{ff2},'pas')
                                    if ~isfield(neuron{n}.mech{t}.range,'pas')  % make range vector if not existing yet
                                        neuron{n}.mech{t}.range.pas.g = NaN(numel(tree{t}.X),1);
                                        neuron{n}.mech{t}.range.pas.cm = NaN(numel(tree{t}.X),1);
                                    end
                                    neuron{n}.mech{t}.range.pas.g(ind) =  neuron{n}.mech{t}.(rangefields{ff1}).pas.g * scale(ind); % add the scaled regional variable only to nodes in that region
                                    neuron{n}.mech{t}.range.pas.cm(ind) =  neuron{n}.mech{t}.(rangefields{ff1}).pas.cm * scale(ind); % add the scaled regional variable only to nodes in that region
                                elseif scaleAll
                                    fields3 = fieldnames(neuron{n}.mech{t}.(rangefields{ff1}).(fields2{ff2}));
                                    for ff3 = 1:numel(fields3)
                                        if ~isempty(strfind(fields3{ff3},'bar'))
                                            if ~isfield(neuron{n}.mech{t}.range,fields2{ff2}) || ~isfield(neuron{n}.mech{t}.range.(fields2{ff2}),fields3{ff3})  % make range vector if not existing yet
                                                neuron{n}.mech{t}.range.(fields2{ff2}).(fields3{ff3}) = NaN(numel(tree{t}.X),1);
                                            end
                                            neuron{n}.mech{t}.range.(fields2{ff2}).(fields3{ff3})(ind) = neuron{n}.mech{t}.(rangefields{ff1}).(fields2{ff2}).(fields3{ff3}) * scale(ind);
                                        end
                                    end
                                end
                            end
                        end
                    else
                        if ~strcmp('all',fields{f1})
                            fields_all = fieldnames(neuron{n}.mech{t}.all);
                            fields_all = fields_all(cellfun(@(x) isempty(strfind(x,'_ion')),fields_all)); % sort out ion mechanisms
                            if ~isempty(setdiff(fields_all,fields2)) % if there are mechanisms defined in "all" regions that have not been addressed in the current region
                                if any(strcmp(fields_all,'spines'))
                                    warning('CAUTION! Field "spines" has been defined in region %s and in region "all"! Please check!',fields{f1})
                                else
                                    morefields = setdiff(fields_all,fields2);
                                    for f2 = 1:numel(morefields)
                                        neuron{n}.mech{t}.(fields{f1}).(morefields{f2}) = neuron{n}.mech{t}.all.(morefields{f2});  % property is now explicitly defined for that region so that it can be modified by spines information
                                    end
                                    fields2 = cat(1,fields2,morefields);
                                end
                                
                            end
                        end
                        for f2 = 1:numel(fields2)
                            if strcmp(fields2{f2},'pas')
                                neuron{n}.mech{t}.(fields{f1}).pas.g = neuron{n}.mech{t}.(fields{f1}).pas.g * scale;
                                neuron{n}.mech{t}.(fields{f1}).pas.cm = neuron{n}.mech{t}.(fields{f1}).pas.cm * scale;
                            elseif scaleAll
                                fields3 = fieldnames(neuron{n}.mech{t}.(fields{f1}).(fields2{f2}));
                                for f3 = 1:numel(fields3)
                                    if ~isempty(strfind(fields3{f3},'bar'))
                                        neuron{n}.mech{t}.(fields{f1}).(fields2{f2}).(fields3{f3}) = neuron{n}.mech{t}.(fields{f1}).(fields2{f2}).(fields3{f3}) * scale;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if isfield(neuron{n},'experiment') && ~isempty(neuron{n}.experiment)
        neuron{n}.experiment = strcat(neuron{n}.experiment,'_scaledspines');
    else
        neuron{n}.experiment = 'Scaledspines';
    end
end

if cellflag
    neuron = neuron{1};
end