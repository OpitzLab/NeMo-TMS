function  [neuron] = t2n_blockchannel(neuron,channels,amount,regions,specify)
% This function manipulates the neuron structure to reduce specific channel
% conductances defined by "channels" in specific regions of the cell(s)
% defined by "regions" by a certain amount defined by "amount"
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms
% channels          string or cell array of strings which channels should
%                   be blocked (i.e. its conductance reduced)
% amount            amount of blockade (0-100) [%]. This could also be used 
%                   to increase conductances if a negative percentage is
%                   given, however it is recommended to rather use
%                   t2n_changemech.
% regions           (optional) allows to select only specific regions of 
%                   cell for blockade. Can be string or cell array of 
%                   strings. The region name has to be defined in the tree 
%                   regions, too.
% specify           (optional) allows to select only specific conductances 
%                   of a channel, eg. gabkbar_BK
%
% OUTPUTS
% neuron            manipulated t2n neuron structure
%
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************



if nargin < 3 || isempty(amount)
    amount = 100;
end
if any(amount > 100)
    error('Values above 100 are not allowed for "amount"')
end

if nargin < 2 || isempty(channels)
    if isempty(channels{1})
        error('No channel to block specified')
    end
elseif ~iscell(channels)
    channels = {channels};
end
if nargin < 5 || isempty(specify)
    specify = {};
elseif ischar(specify) % allows to select only specific conductances of a channel
    specify = {specify};
end
if nargin < 4 || isempty(regions)
    regions = {};
elseif ischar(regions) % allows to select only specific regions of cell
    regions = {regions};
end

if strcmp(channels{1},'except')
    except = 1;
    notchannels = channels;
else
    except = 0;
end
if numel(amount) == 1
    if ~except
        amount = repmat(amount,numel(channels),1);
    end
elseif  numel(amount) ~= numel(channels)
    error('Number of channels to block (m) and number of relative blocking specifications are not the same. Relative blocking variable must have 1 or m elements')
end
flag = false(numel(channels),1);
for t = 1:numel(neuron.mech)
    fields = fieldnames(neuron.mech{t});
    if ~isempty(regions)
        if isfield(fields,'range')
           if ~isempty(intersect(channels,fieldnames(neuron.mech{t}.range)))
                warning('Caution! One or more channels that should be blocked are (at least partly) defined in the "range" variable but block was restricted to a specific region. This could be problematic, if the range variable contains the conductance (e.g. gbar) parameter that should be blocked. This conductance has to be blocked manually then.')
           end
        end
        fields = intersect(fields,regions);
    end
    for f1 = 1:numel(fields)
        fields2 = fieldnames(neuron.mech{t}.(fields{f1}));
        if except
            channels = setdiff(fields2,notchannels);
            amount = repmat(100,numel(channels),1);
        end
        for f2 = 1:numel(fields2)
            for c = 1:numel(channels)
                if ~isempty(strfind(fields2{f2},channels{c}))
                    fields3 = fieldnames(neuron.mech{t}.(fields{f1}).(fields2{f2}));
                    if ~isempty(specify) && numel(specify) >= c && ~isempty(specify{c}) && (ischar(specify{c}) || iscell(specify{c}))  % allows to select only specific conductances of a channel
                        fields3 = intersect(fields3,specify{c});
                        for f3 = 1:numel(fields3)
                            neuron.mech{t}.(fields{f1}).(fields2{f2}).(fields3{f3}) = neuron.mech{t}.(fields{f1}).(fields2{f2}).(fields3{f3}) * (100-amount(c))/100;
                        end
                        flag(c) = true;
                    else
                        for f3 = 1:numel(fields3)
                            if ~isempty(strfind(fields3{f3},'bar')) || (strcmp(channels{c},'pas') && strcmp(fields3{f3},'g'))
                                neuron.mech{t}.(fields{f1}).(fields2{f2}).(fields3{f3}) = neuron.mech{t}.(fields{f1}).(fields2{f2}).(fields3{f3}) * (100-amount(c))/100;
                                flag(c) = true;
                            end
                        end
                    end
                end
            end
        end
    end
end

if ~except && any(~flag)
    warning('Caution,channel(s) %s was/were not found in the mechanism definitions!\n',strcat(channels{~flag}))
end

str = '';
for c = 1:numel(channels)
    str = sprintf('%s%s(%d%%),',str,channels{c},amount(c));
end
str = str(1:end-1);
if isfield(neuron,'experiment') && ~isempty(neuron.experiment)
    neuron.experiment = strcat(neuron.experiment,sprintf('_Block-%s',str));
else
    neuron.experiment = sprintf('Block-%s',str);
end

