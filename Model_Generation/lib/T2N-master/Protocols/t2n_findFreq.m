function [amp] = t2n_findFreq(neuron,tree,desNum,options)
% This function finds the current necessary to let each neuron spike a certain
% amount of spikes with an IClamp protocol previously defined in "neuron" 
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms (see documentation)
% tree              tree cell array with morphologies (see documentation)
% desNum            desired number of spikes during the IClamp protocol
% options           (optional) options for starting T2N
%
% OUTPUT
% amp               current amplitude [nA] for each neuron to reach the desired number of spikes
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 6
    options = '-q-d';
end

counterthresh = 20;


step = 0.05;
neuron.params.tstop = neuron.params.tstop - neuron.pp{1}.IClamp.times(2);
ready = zeros(numel(tree),1);
maxamp = NaN(numel(tree),2) ;
minamp = NaN(numel(tree),2) ;
neuron.params.skiprun = 0;
flag = false;
ostep = zeros(numel(tree),1);
counter = 0;
sflag = ones(numel(tree),1);
if isfield(neuron,'con')
    neuron = rmfield(neuron,'con'); % delete all connections since this is not desired here
end
if isfield(neuron,'record')
    neuron = rmfield(neuron,'record'); % delete all connections since this is not desired here
end
treeind = 1:numel(tree);
while ~flag && counter <= counterthresh
    counter = counter +1;

    for t = 1:numel(tree)
        if ~isfield(tree{t},'artificial')
            %             neuron.record{t}.cell = struct('node',1,'record','v');%;recnode,'ik_Kir';recnode,'gka_ichan3';recnode,'i_pas'};
            if counter == 1
                amp(treeind(t)) = neuron.pp{t}.IClamp.amp(2);
                neuron.pp{t}.IClamp.times = neuron.pp{t}.IClamp.times - neuron.pp{t}.IClamp.times(2)+1;
                neuron.record{t}.cell = struct('node',1,'record','v');
            else
                neuron.pp{t}.IClamp.amp(2) = amp(treeind(t));
            end
            neuron.APCount{t} = [1,-20];
        else
            ready(treeind(t)) = true;
        end
    end
    [out] = t2n(neuron,tree,options);
%     if out.error
%         amp = NaN;
%         return
%     end
    for t = 1:numel(tree)
        if ~isfield(tree{t},'artificial')
            
            thisNum = numel(out.APCtimes{t}{1});
            
            fprintf('Reached: %g spikes of target spikes %g \n',thisNum,desNum)
            if thisNum > desNum  % too many spike
                maxamp(treeind(t),:) = [amp(treeind(t)) thisNum];
                if counter == 1 || sflag(treeind(t))
                    amp(treeind(t)) = amp(treeind(t)) - step;
                else
                    ostep(treeind(t)) = ostep(treeind(t)) / 2;
                    amp(treeind(t)) = amp(treeind(t)) - ostep(treeind(t));
                end
                
                continue
            end
            sflag(treeind(t)) = 0;
            if thisNum == desNum
                ready(treeind(t)) = true;
                if all(ready)
                    flag = true;
                end
            else % lower spikes
                if isnan(maxamp(treeind(t))) || isnan(minamp(treeind(t)))
                    minamp(treeind(t),:) = [amp(treeind(t)) thisNum];
                    ostep(treeind(t)) = sign(desNum-thisNum)*step;
                    amp(treeind(t)) = amp(treeind(t)) + ostep(treeind(t));
                else
                    minamp(treeind(t),:) = [amp(treeind(t)) thisNum];
                    amp(treeind(t)) = (maxamp(treeind(t))-minamp(treeind(t)))/2 + minamp(treeind(t));
                end
            end
        end
    end
    
    [~,~,ib] = intersect(find(ready),treeind);
    %!%!%!
    tree(ib) = [];
    neuron.mech(ib) = [];
    neuron.pp(ib) = [];
    neuron.APCount(ib) = [];
    treeind(ib) = [];
end
if counter > counterthresh
    warndlg('Not all target spikes could be reached! Check!')
end
