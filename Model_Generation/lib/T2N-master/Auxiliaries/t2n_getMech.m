function [varargout] = t2n_getMech(neuron,tree,var)
% This function maps the values of a given NEURON variable onto each
% morphology and optionally returns these values, too;
% INPUTS
% neuron    t2n neuron structure or cell array of structure (see
%           documentation)
% tree      TREES toolbox morphologies
% var       string with valid NEURON variable name, such as 'g_pas'
%
% OUTPUT
% varVec    m-by-n cell array with mapped values for m neuron
%           specifications and n trees
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 3 || isempty(var)
    var = 'g_pas';
    disp('g_pas is used as mechanism')
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
if any(strcmp(var,{'cm','Ra'}))
    var = strcat(var,'_pas');
end
varSplit = regexp(var,'_','split');  % split variable into mechanism name and variable name
varVec = cell(numel(neuron),numel(tree));

stdval = NaN;
if nargin > 3 && exist(fullfile(pwd,'lib_mech',sprintf('%s.mod',mech)),'file')  % check if mod file of mechanism exists
    % read everything
    fid  = fopen(fullfile(pwd,'lib_mech',sprintf('%s.mod',mech)),'r');
    text = textscan(fid,'%s','Delimiter','');
    text = text{1};
    fclose(fid);
    % trim text to PARAMETER part
    LW = regexp(text,'PARAMETER','tokens');
    text = text(find(~cellfun(@isempty,LW)):end);
    LW = regexp(text,'}','tokens');
    text = text(1:find(~cellfun(@isempty,LW),1,'first'));
    % find lines that match the par name
    LW = regexp(text,par,'tokens');
    if sum(~cellfun(@isempty,LW)) > 1
        warning('Parameter name occurs more than once in the PARAMETER section of the mod file. Hence, standard parameter value cannot be read. Delete commented lines for proper reading of value.')
    else
        text = text{~cellfun(@isempty,LW)};
        LW = regexp(text,[par,'[\s\.=]+[+-]?(\d+)[\.]?(\d+)[\r\t\s]'],'tokens');  % extract numeric value
        LW = sprintf('%s.',LW{1}{:});  % stitch extracted number strings with decimal point
        stdval = str2double(LW(1:end-1)); % make number
    end
end

for n = 1:numel(neuron) % go through all neuron definitions
    for t = 1:numel(tree)  % go through all morphologies
        varVec{n,t} = ones(numel(tree{t}.X),1)*stdval;  % initialize the vector as nans or the standard value from the mod file
        if isfield(neuron{n},'mech')   % check for mechanism definition
            regions = fieldnames(neuron{n}.mech{t}); % get region names
            if any(strcmp(regions,'all'))   % variable was set in all nodes
                if isfield(neuron{n}.mech{t}.all,varSplit{2})  % check if mechanism exists
                    if isfield(neuron{n}.mech{t}.all.(varSplit{2}),varSplit{1}) % check if variable is defined at that location
                        varVec{n,t}(:) = neuron{n}.mech{t}.all.(varSplit{2}).(varSplit{1}); % save value in all nodes
                    elseif isnan(stdval)   % check if standard value was not found
                        varVec{n,t}(:) = Inf;  % mechanism seems to be defined but not the variable. remember that
                    end
                end
               regions = regions(~strcmp(regions,'all')); % delete region 'all' from regions
            end
            regions = intersect(regions,tree{t}.rnames);  % remove all nondefined regions
            for r = 1:numel(regions)  % go through regions
                if isfield(neuron{n}.mech{t}.(regions{r}),varSplit{2})  % check if mechanism exists at that region
                    ind = tree{t}.R == find(strcmp(tree{t}.rnames,regions{r}));  % get index to all nodes that belong to this region
                    if isfield(neuron{n}.mech{t}.(regions{r}).(varSplit{2}),varSplit{1}) % check if variable is defined at that location
                        varVec{n,t}(ind) = neuron{n}.mech{t}.(regions{r}).(varSplit{2}).(varSplit{1}); % save value in these nodes
                    elseif isnan(stdval)  % check if standard value was not found
                        varVec{n,t}(ind) = Inf;  % mechanism seems to be defined but not the variable. remember that
                    end
                end
            end
            if isfield(neuron{n}.mech{t},'range')  % check for a range structure
                mechanisms = intersect(varSplit{2},fieldnames(neuron{n}.mech{t}.range)); % check if our mechanism exists in the range variable
                if ~isempty(mechanisms)
                    vars = intersect(varSplit{1},fieldnames(neuron{n}.mech{t}.range.(mechanisms{1}))); % check if our variable exists in the range variable
                    if ~isempty(vars)
                        ind = ~isnan(neuron{n}.mech{t}.range.(mechanisms{1}).(vars{1}));  % get only the not nan values (nan means use the default value or value defined by the region spec)
                        varVec{n,t}(ind) = neuron{n}.mech{t}.range.(mechanisms{1}).(vars{1})(ind); % save the not nan values
                    end
                end
            end
        end
        if any(isinf(varVec{n,t}))
            warning('Caution! Values of %s seem to have not been defined at some nodes! NEURON uses the mechanisms default value.',var)
            varVec{n,t}(isinf(varVec{n,t})) = NaN;
        end
        if nargout == 0  % map values on tree
            figure
            if isfield(tree{t},'name')
                tname = tree{t}.name;
            else
                tname = sprintf('tree %d',t);
            end
            title(sprintf('Neuron Simulation: %d, tree: "%s", variable "%s"',n,strrep(tname,'_','\_'),strrep(var,'_','\_')))
            plot_tree(tree{t},varVec{n,t});
            lims = [min(varVec{n,t}),max(varVec{n,t})];
            if diff(lims) == 0  % limits HAVE to be different from each other
                lims = [lims(1)-1e-9,lims(1)+1e-9];
            end
            if ~all(isnan(lims))
                set(gca,'CLim',lims)
            end
        end
    end
end
if nargout > 0
    varargout{1} = varVec;
end
end

