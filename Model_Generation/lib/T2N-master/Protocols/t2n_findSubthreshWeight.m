function [weight] = t2n_findSubthreshWeight(neuron,tree,weight,freq,tim)
% This function finds the synaptic weight of Exp2Syn synapses necessary to 
% have the neurons defined in "neuron" and "tree" at a sub-spiking threshold level
% 
% INPUT
% neuron	T2N neuron structure (already containing all synapses for which
%           the weight should be searched
% tree      TREES toolbox tree cell array
% weight	starting weight for the synapses
% freq      frequency of the netstim input
% tim (optional): time [ms] after which simulation stops (default five
% rounds of stimulation)
%
% OUTPUT
% weight	synaptic weight for each neuron to reach subthreshold state
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 5 || isempty(weight)
    weight = ones(numel(tree),1) * 0.001;
elseif numel(weight)==1
    weight = ones(numel(tree),1) * weight;
end
if nargin < 6 || isempty(freq)
    freq = 10; % Hz
end
if nargin < 7 || isempty(tim)
    tim = 1/freq*1000*5; % 5 times input 1/frequency
end
counterthresh = 20;


e = 0.00005;
ostep = ones(numel(tree),1)*0.0001;

neuron.params.cvode = 0;
neuron.params.dt = 0.5; 
neuron.params.tstop = tim + 10;

ready = zeros(numel(tree),1);
spikeamp = NaN(numel(tree),1) ;
neuron.params.skiprun = 0;
flag = false;

counter = 0;
sflag = true(numel(tree),1);
if iscell(neuron)
    neuron = neuron{1};
end
if isfield(neuron,'con')
    neuron = rmfield(neuron,'con'); % delete all connections since this is not desired here
end
tree{end+1} = struct('artificial','NetStim','start',10,'interval',1/freq*1000,'number',10);
        
treeind = 1:numel(tree);
while ~flag && counter <= counterthresh
    counter = counter +1;

    for t = 1:numel(tree)
        if ~isfield(tree{t},'artificial')
            neuron.record{t}.cell = struct('node',1,'record','v');
            
            neuron.con(t) = struct('source',struct('cell',numel(tree),'watch','on'),'target',struct('cell',t,'pp','Exp2Syn','node',neuron.pp{t}.Exp2Syn.node),'weight',weight(treeind(t)),'delay',0,'threshold',0.5);
            neuron.APCount{t} = [1,-30];
        else
            ready(treeind(t)) = 2;
        end
    end
    [out] = t2n(neuron,tree,'-q-d');
    if out.error
        weight = NaN;
        return
    end
    for t = 1:numel(tree)
        if ~isfield(tree{t},'artificial')
                      
            maxv = max(out.record{t}.cell.v{1});
            if maxv > -40 || ~isempty(out.APCtimes{t}{1})  % spike
                spikeamp(treeind(t),:) = weight(treeind(t));
                if abs(ostep(treeind(t))) < e % numel(out.APCtimes{t}{1}) == 1
                    ready(treeind(t)) = true;
                    if all(ready)
                        flag = true;
                    end
                end
                fprintf('Reached Spike with %g microS weight\n',weight(treeind(t)));
                if counter == 1 || sflag(treeind(t))
                    ostep(treeind(t)) = min(weight(treeind(t)),ostep(treeind(t)) * 2);
                    weight(treeind(t)) = weight(treeind(t)) - ostep(treeind(t));
                else
                    ostep(treeind(t)) = ostep(treeind(t)) / 2;
                    weight(treeind(t)) = weight(treeind(t)) - ostep(treeind(t));
                end
                
                continue
                
            end
            sflag(treeind(t)) = 0;
            fprintf('Reached %g mV with %g microS weight\n',maxv,weight(treeind(t)));

            if ~isnan(spikeamp(treeind(t))) && abs(spikeamp(treeind(t)) - weight(treeind(t))) < e
                weight(treeind(t)) = spikeamp(treeind(t));
                ready(treeind(t)) = true;
                if all(ready)
                    flag = true;
                end
            else
                if isnan(spikeamp(treeind(t)))
                    weight(treeind(t)) = weight(treeind(t)) + ostep(treeind(t));
                else
                    ostep(treeind(t)) = ostep(treeind(t)) / 2;
                    weight(treeind(t)) = weight(treeind(t)) + ostep(treeind(t));
                end
            end
        end
    end
    
    [~,~,ib] = intersect(find(ready==1),treeind);
    
    if ~isempty(ib)
        neuron.con(ib) = [];
        tree(ib) = [];
        neuron.mech(ib) = [];
        neuron.record(ib) = [];
        neuron.pp(ib) = [];
        neuron.APCount(ib) = [];
        treeind(ib) = [];
    end
end
if counter > counterthresh
    warndlg('Not all target voltages could be reached! Check!')
end
