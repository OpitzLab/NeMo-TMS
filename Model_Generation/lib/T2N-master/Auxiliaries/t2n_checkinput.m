function [neuron,tree,usestreesof,nocell,exchfolder, tag2Node, legacyOutput] = t2n_checkinput(neuron,tree)
% This function checks the neuron structure for correct definition of the
% used morphologies and returns info about it
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms (see documentation)
% tree              tree cell array with morphologies (see documentation)
%
% OUTPUTS
% tree              corrected tree cell array
% neuron            corrected neuron structure
% usestreesof       points to the neuron entry/instance from which the tree
%                   definitions are taken from
% nocell            Boolean if neuron input was a structure or cell array
% exchfolder        name for exchfolder that was possibly found in the
%                   neuron structure
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if ~exist(fullfile(pwd,'morphos','hocs'),'dir')
    mkdir(fullfile(pwd,'morphos','hocs'));
end

%% check neuron
nocell = false;
if isstruct(neuron)         % transform structure neuron to cell neuron
    if numel(neuron) == 1
        neuron = {neuron};
        nocell = true;
    else
        neuron = arrayfun(@(y) {y},neuron);
    end
end

%% check tree
if iscell(tree) && iscell(tree{1})
    tree = tree{1};
elseif isstruct(tree)
    tree = {tree};
end
for t = 1:numel(tree)
    if isfield(tree{t},'artificial') && ~isfield(tree{t},'NID')
        tree{t}.NID = strcat('cell_',tree{t}.artificial);           % artificial cells only need one morph hoc file which is named cell_ + the name of the artificial cell..
    end
end
if  ~all(cellfun(@(x) isfield(x,'NID'),tree))
    doit = 1;
else
    NIDs = unique(cellfun(@(x) x.NID,tree,'UniformOutput',0));  % improves speed if many same cells are used
    doit = 0;
end
if doit || ~all(cellfun(@(x) exist(fullfile(pwd,'morphos','hocs',strcat(x,'.hoc')),'file'),NIDs))
    answer = questdlg('Caution! Not all of your trees have been transformed for NEURON yet or hoc file is missing! Transforming now..','Transform trees','OK','Cancel','OK');
    if strcmp(answer,'OK')
        ind = cellfun(@(x) ~isfield(x,'NID'),tree) ;
        if ~all(ind)
            ind(~ind) = ~cellfun(@(x) exist(fullfile(pwd,'morphos','hocs',strcat(x.NID,'.hoc')),'file'),tree(~ind));
        end
        tree(ind) = t2n_writeTrees(tree(ind));
    else
        error('Aborted');
    end
end

%% check tree/neuron consistency
thesetrees = cell(numel(neuron),1);
usestreesof = zeros(numel(neuron),1);
flag = false;
bool = cellfun(@(y) isfield(y,'tree'),neuron);
if all(cellfun(@(y) ischar(y.tree),neuron(bool))) % no trees defined (only due to t2n_as function). make bool = 0
    bool(:) = 0;
end
nempty = cellfun(@isempty,neuron);
if any(nempty)
    error('The defined simulation #%d is empty, please check\n',find(nempty))
end
switch sum(bool)
    case 1   % use that treeids defined in that one simulation
        if isnumeric(neuron{bool}.tree)
            thesetrees = repmat({unique(neuron{bool}.tree)},numel(neuron),1);
            usestreesof = repmat(find(bool),numel(neuron),1);
        else
            x = t2n_getref(1,neuron,'tree');
            if isempty(x) % means it refers to itself (maybe due to usage of t2n_as)..use normal trees..
                thesetrees = repmat({1:numel(tree)},numel(neuron),1);
                usestreesof = ones(numel(neuron),1);
            else
                n = find(bool);
                flag = true;
            end
        end
    case 0      % if no trees are given, use trees that are given to t2n in their order...
        thesetrees = repmat({1:numel(tree)},numel(neuron),1);
        usestreesof = ones(numel(neuron),1);
    case numel(neuron)
        for n = 1:numel(neuron)
            x = t2n_getref(n,neuron,'tree');
            if ~isnan(x)
                thesetrees{n} = unique(neuron{x}.tree);
                usestreesof(n) = x;
            elseif isempty(x)
                thesetrees{n} = 1:numel(tree);%repmat({1:numel(tree)},numel(neuron),1);
                usestreesof(n) = 1;%ones(numel(neuron),1);
            else
                flag = true;
                break
            end
        end
    otherwise  % if more than one are given, t2n cannot know which trees you want
        n = find(bool);
        flag = true;
end
if flag
    error('Error in neuron{%d}.tree, please check\n',n)
end
refP = t2n_getref(1,neuron,'params');
if ~isnan(refP) && isfield(neuron{refP}.params,'exchfolder')
    exchfolder = neuron{1}.params.exchfolder;
else
    exchfolder = [];
end
tag2Node = cell(numel(neuron),1);
legacyOutput = false(numel(neuron),1);
for n = 1:numel(neuron)
    neuron{n}.tree = thesetrees{n};
    
    %% check for several standard parameter and initialize default value if not set
    refP = t2n_getref(n,neuron,'params');
    refPP = t2n_getref(n,neuron,'pp');
    refR = t2n_getref(n,neuron,'record');
    refPL = t2n_getref(n,neuron,'play');
    if n == refP % only check if current instance has its own parameter struct
        if ~isfield(neuron{n}.params,'parallel')
            neuron{n}.params.parallel = 0;
        end
        if neuron{n}.params.parallel == 1
            warning('Neuron instance %d has params.parallel set to 1, however this variable is no boolean but defines the number of cores that should be used. Using 1 core has not advantage.',n)
        elseif numel(neuron) > 1 && neuron{n}.params.parallel > 0
            warning('It seems you have enabled parallel Neuron together with multiple Neuron instances. Be sure that there are %d cores available, otherwise there will be no improvement in speed and the CPUs might be overloaded.',numel(neuron)*neuron{n}.params.parallel)
        end
        if ~isfield(neuron{n}.params,'cvode')
            neuron{n}.params.cvode = false;
        end
        if ~isfield(neuron{n}.params,'use_local_dt')
            neuron{n}.params.use_local_dt = 0;
        end
        if ~isfield(neuron{n}.params,'nseg')
            neuron{n}.params.nseg = 'dlambda';
            disp('Number of segments or nseg rule not set in neuron{n}.params.nseg. Dlambda rule will be applied')
        elseif strcmpi(neuron{n}.params.nseg,'d_lambda')
            neuron{n}.params.nseg = 'dlambda';
        end
        if ~isfield(neuron{n}.params,'dlambda')
            neuron{n}.params.dlambda = 0.1;
        end
        if ~isfield(neuron{n}.params,'freq')
            neuron{n}.params.freq = 300;
        end
        if ~isfield(neuron{n}.params,'tstart')
            neuron{n}.params.tstart = 0;
        end
        if ~isfield(neuron{n}.params,'tstop')
            neuron{n}.params.tstop = 200;
            disp('Tstop not defined in neuron{n}.params.tstop. Default value of 200 ms is applied.')
        end
        if ~isfield(neuron{n}.params,'dt')
            neuron{n}.params.dt = 0.025;
            if ~neuron{n}.params.cvode
                disp('Time step not defined in neuron{n}.params.dt. Default value of 0.025 ms is applied.')
            end
        end
        if ~isfield(neuron{n}.params,'accuracy')
            neuron{n}.params.accuracy = 0;
        end
        if ~isfield(neuron{n}.params,'skiprun')
            neuron{n}.params.skiprun = false;
        end
        if ~isfield(neuron{n}.params,'q10')
            neuron{n}.params.q10 = false;
        end
        if ~isfield(neuron{n}.params,'prerun')
            neuron{n}.params.prerun = false;
        end
        if neuron{n}.params.cvode && isnumeric(neuron{n}.params.dt) && (~isnan(neuron{n}.params.dt) || ~isempty(neuron{n}.params.dt))
            warning ('t2n:cvode', 'Dt is set but cvode is active. Dt will be ignored');
        end
    end
    % stupid workaround if someone put the node information along second
    % dimension...
    if refPP == n
        allTags = [];
        allNodes = [];
        for p = 1:numel(neuron{n}.pp)
            fields = fieldnames(neuron{n}.pp{p});
            for f = 1:numel(fields)
                for nn = 1:numel(neuron{n}.pp{p}.(fields{f}))
                    neuron{n}.pp{p}.(fields{f})(nn).node = neuron{n}.pp{p}.(fields{f})(nn).node(:);
                    if ~isfield(neuron{n}.pp{p}.(fields{f})(nn),'tag') || isempty(neuron{n}.pp{p}.(fields{f})(nn).tag)
                        % create a unique tag for those who haven't any
                        neuron{n}.pp{p}.(fields{f})(nn).tag = string(arrayfun(@(x) sprintf('tag_%d_%s_%d_%d_%d',p,fields{f},nn,neuron{n}.pp{p}.(fields{f})(nn).node(x),x),1:numel(neuron{n}.pp{p}.(fields{f})(nn).node),'uni',0));
                    elseif numel(neuron{n}.pp{p}.(fields{f})(nn).node) ~= numel(neuron{n}.pp{p}.(fields{f})(nn).tag)
                        error('t2n:InconsistentTag','The number of specified tags does not equal the number of nodes specified for instance %d of %s (cell %d, neuron instance %d). It might be possible that you specified tags in a cell array of strings. In this case please convert to string array with string(X) before defining the point process!',nn,fields{f},p,n)
                    end
                    allTags = cat(1,allTags,neuron{n}.pp{p}.(fields{f})(nn).tag);
                    allNodes = cat(1,allNodes,neuron{n}.pp{p}.(fields{f})(nn).node);
                end
            end
        end
        if numel(allTags) ~= numel(unique(allTags))
            error('t2n:TagDuplicates','The point process tag strings specified have to be unique in neuron instance %d! Maybe you added multiple point processes of the same type to the same node without applying a unique tag?',n)
        end
        tag2Node{n} = cell2struct(num2cell(allNodes), cellstr(allTags));
    end
    
    % check recording structure
    if refR == n
        for tt = 1:numel(neuron{n}.tree)
            t = neuron{n}.tree(tt);
            if numel(neuron{refR}.record) >= t && ~isempty(neuron{refR}.record{t})  % if a recording site was defined for  this tree
                recfields = fieldnames(neuron{refR}.record{t});
                if isfield(tree{neuron{n}.tree(tt)},'artificial')
                    if numel(recfields) > 1 && any(strcmp(recfields,'record'))
                        % put the structure in field named as the artificial neuron
                        neuron{refR}.record{t} = struct(tree{neuron{n}.tree(tt)}.artificial,neuron{refR}.record{t});
                    end
                    recfields = fieldnames(neuron{refR}.record{t});
                end
                
                for f1 = 1:numel(recfields)
                    % these lines filter out multiple
                    % recording node definitions
                    if numel(setdiff(fieldnames(neuron{refR}.record{t}.(recfields{f1})),{'tag','record','node','tvec','Dt'}))>0
                        error('this has not been implemented yet. write to marcel.beining@gmail.com with specification of the error line')
                    end
                    % find recording fields with same recording variable
                    [uniqrecs,~,indrecgroups] = unique({neuron{refR}.record{t}.(recfields{f1}).record});
                    % rearrange rec structure such that each recording
                    % is in its own structure. Leave this out if
                    % several groups of pps are defined
                    if strcmp(recfields{f1},'cell') || numel(neuron{refPP}.pp{t}.(recfields{f1})) == 1
                        tmpstruct = neuron{refR}.record{t}.(recfields{f1})([]);
                        for u = 1:numel(uniqrecs)  % go through variable groups
                            if isfield(neuron{refR}.record{t}.(recfields{f1}),'node')
                                % get the unique nodes for that recorded variable
                                unodes = unique(cat(1,neuron{refR}.record{t}.(recfields{f1})(indrecgroups==u).node));
                                % save these in a temporary structure
                                tmpstruct(u).node = unodes;
                            else
                                % get the unique tags for that recorded variable
                                if numel(neuron{refR}.record{t}.(recfields{f1})) > 1
                                    utags = unique(string(cat(1,neuron{refR}.record{t}.(recfields{f1})(indrecgroups==u).tag)));
                                else
                                    utags = string(neuron{refR}.record{t}.(recfields{f1}).tag);
                                end
                                % save these in a temporary structure
                                tmpstruct(u).tag = utags;
                            end
                            tmpstruct(u).record = uniqrecs{u};
                            
                            if isfield(neuron{refR}.record{t}.(recfields{f1}),'Dt')
                                tmpstruct(u).Dt = unique(cat(1,neuron{refR}.record{t}.(recfields{f1})(indrecgroups==u).Dt));
                                if numel(tmpstruct(u).Dt) > 1
                                    error('Different resampling values Dt given for the same recording variable %s. This is not supported!',uniqrecs{u})
                                end
                            end
                            if isfield(neuron{refR}.record{t}.(recfields{f1}),'tvec')
                                tmpstruct(u).tvec = unique(cat(1,neuron{refR}.record{t}.(recfields{f1})(indrecgroups==u).tvec));
                                if numel(tmpstruct(u).tvec) > 1
                                    error('Different resampling vectors tvec given for the same recording variable %s. This is not supported!',uniqrecs{u})
                                end
                            end
                        end
                        neuron{refR}.record{t}.(recfields{f1}) = tmpstruct;  % overwrite old record defiition with new record structure
                    else
                        warning('It seems recordings of different PP groups have been defined. Make sure that indices match, e.g. .record{1}.ExpSyn(3) is to target only .pp{1}.ExpSyn(3) etc.')
                    end
                    
                    % check if recording is a pp, then replace node
                    % indices with tag indices
                    if ~isfield(tree{neuron{n}.tree(tt)},'artificial') && ~strcmp(recfields{f1},'cell')
                        for r = 1:numel(neuron{refR}.record{t}.(recfields{f1}))
                            if isfield(neuron{refR}.record{t}.(recfields{f1}),'node')
                                if isfield(neuron{refR}.record{t}.(recfields{f1})(r),'tag') && ~isempty(neuron{refR}.record{t}.(recfields{f1})(r).tag)
                                    warning('You defined tags and nodes in record for %s. Node entry is removed',recfields{f1})
                                else
                                    % check if ppg feature was used
                                    if isfield(neuron{refR}.record{t}.(recfields{f1})(r),'ppg')
                                        ppg = neuron{refR}.record{t}.(recfields{f1})(r).ppg;
                                    elseif numel(neuron{refPP}.pp{t}.(recfields{f1})) == numel(neuron{refR}.record{t}.(recfields{f1}))
                                        ppg = r;
                                    else
                                        ppg =1;
                                    end
                                    tags = cell(numel(neuron{refR}.record{t}.(recfields{f1})(r).node),1);
                                    % find pp tags that correspond to the record node
                                    % indices
                                    for in =  1:numel(neuron{refR}.record{t}.(recfields{f1})(r).node)
                                        ind = find(neuron{refPP}.pp{t}.(recfields{f1})(ppg).node == neuron{refR}.record{t}.(recfields{f1})(r).node(in));
                                        if numel(ind) > 1
                                            warning('There is more than one %s in tree %d at node %d but you did not use the ppg or tag feature to distinguish. Only from one %s is recorded now!',recfields{f1},t,neuron{refR}.record{t}.(recfields{f1})(r).node(in),recfields{f1})
                                            ind = ind(1);
                                        elseif isempty(ind)
                                            warning('Node %d of cell %d does not comprise the PP "%s". Recording is ignored.',neuron{refR}.record{t}.(recfields{f1})(r).node(in),t,recfields{f1})
                                            continue
                                        end
                                        tags{in} = neuron{refPP}.pp{t}.(recfields{f1})(r).tag(ind);
                                    end
                                    tags = tags(~cellfun(@isempty,tags));
                                    neuron{refR}.record{t}.(recfields{f1})(r).tag = string(tags);  
                                    legacyOutput(n) = true;
                                end

                            end
                        end
                        % remove node field
                        neuron{refR}.record{t}.(recfields{f1}) = rmfield(neuron{refR}.record{t}.(recfields{f1}),'node');
                    end
                end
            end
        end
    else
        legacyOutput(n) = legacyOutput(refR);
    end
    
    % check play structure
    if refPL == n
        for tt = 1:numel(neuron{n}.tree)
            t = neuron{n}.tree(tt);
            if numel(neuron{refPL}.play) >= t && ~isempty(neuron{refPL}.play{t})  % if a recording site was defined for  this tree
                playfields = fieldnames(neuron{refPL}.play{t});
                if isfield(tree{neuron{n}.tree(tt)},'artificial')
                    if numel(playfields) > 1 && any(strcmp(playfields,'play'))
                        % put the structure in field named as the artificial neuron
                        neuron{refPL}.play{t} = struct(tree{neuron{n}.tree(tt)}.artificial,neuron{refPL}.play{t});
                    end
                end
                playfields = fieldnames(neuron{refPL}.play{t});
                
                for f1 = 1:numel(playfields)
                    % these lines filter out multiple
                    % playing node definitions
                    if numel(setdiff(fieldnames(neuron{refPL}.play{t}.(playfields{f1})),{'tag','play','node','times','value','cont'}))>0
                        error('this play field has not been implemented yet. write to marcel.beining@gmail.com with specification of the error line')
                    end
                    % find recording fields with same recording variable
                    [uniqplays,~,indplaygroups] = unique({neuron{refPL}.play{t}.(playfields{f1}).play});
                    % rearrange rec structure such that each recording
                    % is in its own structure. Leave this out if
                    % several groups of pps are defined
                    if strcmp(playfields{f1},'cell') || numel(neuron{refPP}.pp{t}.(playfields{f1})) == 1
                        tmpstruct = neuron{refPL}.play{t}.(playfields{f1})([]);
                        addFields = setdiff(fieldnames(neuron{refPL}.play{t}.(playfields{f1})),["node","tag","play"]);
                        for u = 1:numel(uniqplays)  % go through variable groups
                            if isfield(neuron{refPL}.play{t}.(playfields{f1}),'node')
                                % get the unique nodes for that recorded variable
                                unodes = unique(cat(1,neuron{refPL}.play{t}.(playfields{f1})(indplaygroups==u).node));
                                % save these in a temporary structure
                                tmpstruct(u).node = unodes;
                            else
                                % get the unique tags for that recorded variable
                                if numel(neuron{refPL}.play{t}.(playfields{f1})) > 1
                                    utags = unique(string(cat(1,neuron{refPL}.play{t}.(playfields{f1})(indplaygroups==u).tag)));
                                else
                                    utags = string(neuron{refPL}.play{t}.(playfields{f1}).tag);
                                end
                                % save these in a temporary structure
                                tmpstruct(u).tag = utags;
                            end
                            tmpstruct(u).play = uniqplays{u};
                            for a = 1:numel(addFields)
                                tmpstruct(u).(addFields(a)) = neuron{refPL}.play{t}.(playfields{f1}).(addFields(a));
                            end
                        end
                        neuron{refPL}.play{t}.(playfields{f1}) = tmpstruct;  % overwrite old play definition with new play structure
                    else
                        warning('It seems plays into different PP groups have been defined. Make sure that indices match, e.g. .play{1}.ExpSyn(3) is to target only .pp{1}.ExpSyn(3) etc.')
                    end
                    
                    % check if play target is a pp, then replace node
                    % indices with tag indices
                    if ~strcmp(playfields{f1},'cell')
                        for r = 1:numel(neuron{refPL}.play{t}.(playfields{f1}))
                            if isfield(neuron{refPL}.play{t}.(playfields{f1}),'node')
                                if isfield(neuron{refPL}.play{t}.(playfields{f1}),'tag')
                                    warning('You defined tags and nodes in play for %s. Node entry is removed',playfields{f1})
                                else
                                    % check if ppg feature was used
                                    if isfield(neuron{refPL}.play{t}.(playfields{f1})(r),'ppg')
                                        ppg = neuron{refPL}.play{t}.(playfields{f1})(r).ppg;
                                    elseif numel(neuron{refPP}.pp{t}.(playfields{f1})) == numel(neuron{refPL}.play{t}.(playfields{f1}))
                                        ppg = r;
                                    else
                                        ppg =1;
                                    end
                                    tags = cell(numel(neuron{refPL}.play{t}.(playfields{f1})(r).node),1);
                                    % find pp tags that correspond to the play node
                                    % indices
                                    for in =  1:numel(neuron{refPL}.play{t}.(playfields{f1})(r).node)
                                        ind = find(neuron{refPP}.pp{t}.(playfields{f1})(ppg).node == neuron{refPL}.play{t}.(playfields{f1})(r).node(in));
                                        if numel(ind) > 1
                                            warning('There is more than one %s in tree %d at node %d but you did not use the ppg or tag feature to distinguish. Only to one %s is played now!',playfields{f1},t,neuron{refPL}.play{t}.(playfields{f1})(r).node(in),playfields{f1})
                                            ind = ind(1);
                                        elseif isempty(ind)
                                            warning('Node %d of cell %d does not comprise the PP "%s". Playing is ignored.',neuron{refPL}.play{t}.(playfields{f1})(r).node(in),t,playfields{f1})
                                            continue
                                        end
                                        tags{in} = neuron{refPP}.pp{t}.(playfields{f1})(r).tag(ind);
                                    end
                                    tags = tags(~cellfun(@isempty,tags));
                                    neuron{refPL}.play{t}.(playfields{f1})(r).tag = string(tags);  
                                end
                                % remove node field
                                neuron{refPL}.play{t}.(playfields{f1}) = rmfield(neuron{refPL}.play{t}.(playfields{f1}),'node');
                            end
                        end
                        
                    end
                end
            end
        end
    end
    
    % check con structure
    if t2n_getref(n,neuron,'con') == n
        % check for all con fields and provide standard values if not
        % existent
        for c = 1:numel(neuron{n}.con)
            if ~isfield(neuron{n}.con(c).source,'watch')
                if isfield(tree{neuron{n}.tree==neuron{n}.con(c).source.cell},'artificial') || (isfield(neuron{n}.con(c).source,'pp') && ~isempty(neuron{n}.con(c).source.pp))
                    neuron{n}.con(c).source.watch = 'on';
                else
                    neuron{n}.con(c).source.watch = 'v';
                end
            end
            if neuron{n}.params.parallel && isfield(neuron{n}.con(c).source,'pp')
                error('Defining a point process as the source of a NetCon is not implemented in t2n when parallel mode is activated. Please contact the developer')
                % what would be needed: change reg_cell to accept pps,
                % assure that pps gid is also registered on the same node
                % as its host cell
            elseif ~isfield(neuron{n}.con(c).source,'pp')
                neuron{n}.con(c).source.pp = [];
            end
            if ~isfield(neuron{n}.con(c).source,'node')
                if isfield(tree{neuron{n}.tree==neuron{n}.con(c).source.cell},'artificial')
                    neuron{n}.con(c).source.node = [];
                else
                    neuron{n}.con(c).source.node = 1;
                end
            elseif isfield(neuron{n}.con(c).source,'tag')
                warning('Tags and nodes specified at the same time in NetCons. Only tags are considered')
            end
            if ~isfield(neuron{n}.con(c).source,'ppg')
                if ~isempty(neuron{n}.con(c).source.pp)
                    neuron{n}.con(c).source.ppg = find(arrayfun(@(x) any(ismember(neuron{n}.con(c).source.node,x.node)), neuron{refPP}.pp{neuron{n}.con(c).source.cell}.(neuron{n}.con(c).source.pp))); % connect to all instances of the point process which exist at that source node
                else
                    neuron{n}.con(c).source.ppg = 0;
                end
            end
        end
        % if multiple nodes/cells/pps have been defined at once in the con
        % list, make them single
        if isfield(neuron{n}.con(1).source,'tag')
            neuron{n}.con = detangleCon(neuron{n}.con,'source','tag');
        else
            neuron{n}.con = detangleCon(neuron{n}.con,'source','cell');
            neuron{n}.con = detangleCon(neuron{n}.con,'source','node');
            if ~neuron{n}.params.parallel                % change this is if implemented
                neuron{n}.con = detangleCon(neuron{n}.con,'source','pp');
                neuron{n}.con = detangleCon(neuron{n}.con,'source','ppg');
            end
        end
        if isfield(neuron{n}.con(1).target,'tag')
            neuron{n}.con = detangleCon(neuron{n}.con,'target','tag');
        else
            neuron{n}.con = detangleCon(neuron{n}.con,'target','cell');
            neuron{n}.con = detangleCon(neuron{n}.con,'target','node');
            neuron{n}.con = detangleCon(neuron{n}.con,'target','pp');
        end
    end
end
end

function newcon = detangleCon(oldcon,field1,field2,field3)

newcon = cell(numel(oldcon),1);
for c = 1:numel(oldcon)
    if isfield(oldcon(c).(field1),field2) && ~ischar(oldcon(c).(field1).(field2)) && ~isempty(oldcon(c).(field1).(field2))
        num = numel(oldcon(c).(field1).(field2));
        newcon{c} = repmat(oldcon(c),num,1);
        for m = 1:num
            newcon{c}(m).(field1).(field2) = newcon{c}(m).(field1).(field2)(m);
            if exist('field3','var')
                newcon{c}(m).(field1).(field3) = newcon{c}(m).(field1).(field3)(m);
            end
        end
    else
        newcon{c} = oldcon(c);
    end
end
newcon = cat(1,newcon{:});

end
