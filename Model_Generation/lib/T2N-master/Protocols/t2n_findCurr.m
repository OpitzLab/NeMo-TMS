function [amp, Vrest] = t2n_findCurr(neuron,tree,desv,amp,options)
% This function finds the current necessary to keep the neurons defined in
% "neuron" and "tree" at a desired voltage or alternatively finds the
% spiking threshold.
% 
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms (see documentation)
% tree              tree cell array with morphologies (see documentation)
% desv              desired voltage (mV) or 'spike' if spiking threshold is searched
% amp               (optional) starting values of the current amplitudes (one for each
%                   neuron defined in tree)
% options           (optional) options for starting T2N (see t2n)
%
% OUTPUTS
% amp               current amplitude for each neuron to reach the desired voltage/spike
% Vrest             resting potential of each cell
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 6
    options = '-q-d';
end
if nargin < 5 || isempty(amp)
    amp = zeros(numel(tree),1);
end
counterthresh = 20;

if ischar(desv) && strcmpi(desv,'Spike')
    e = 0.0005; %nA allowed difference to find spiking thresh
    neuron.params.dt = 1;
    neuron.params.tstop = 500;
else
    e = 0.2; % mV allowed difference between target voltage and actual voltage
    neuron.params.dt = 10;
    neuron.params.tstop = 2000;
end
step = 0.05;

neuron.params.cvode = 0;

ready = zeros(numel(tree),1);
maxamp = NaN(numel(tree),2) ;
minamp = NaN(numel(tree),2) ;

flag = false;
ostep = zeros(numel(tree),1);
counter = 0;
Vrest = NaN(numel(tree),1);
sflag = ones(numel(tree),1);
if isfield(neuron,'con')
    neuron = rmfield(neuron,'con'); % delete all connections since this is not desired here
end
if isfield(neuron,'play')
    neuron = rmfield(neuron,'play'); % delete all play vectors since this is not desired here
end
if isfield(neuron,'pp')
    neuron = rmfield(neuron,'pp'); % delete all point processes since this is not desired here
end

treeind = 1:numel(tree);
while ~flag && counter <= counterthresh
    counter = counter +1;

    for t = 1:numel(tree)
        if ~isfield(tree{t},'artificial')
            neuron.record{t}.cell = struct('node',1,'record','v');%;recnode,'ik_Kir';recnode,'gka_ichan3';recnode,'i_pas'};
            neuron.pp{t}.IClamp = struct('node',1,'del',0,'dur',2000,'amp', amp(treeind(t))); %n,del,dur,amp
            neuron.APCount{t} = [1,-20];
        else
            ready(treeind(t)) = true;
        end
    end
    [out] = t2n(neuron,tree,options);

    for t = 1:numel(tree)
        if ~isfield(tree{t},'artificial')
            
            thisv = out.record{t}.cell.v{1}(end);
            
            maxv = max(out.record{t}.cell.v{1});
            if maxv > -40 || ~isempty(out.APCtimes{t}{1})  % spike
                maxamp(treeind(t),:) = [amp(treeind(t)) maxv];
                if ischar(desv) && strcmpi(desv,'Spike')
                    if isnan(minamp(treeind(t)))
                        ostep(treeind(t)) = -step;
                        amp(treeind(t)) = amp(treeind(t)) + ostep(treeind(t));
                    elseif abs(minamp(treeind(t),1)-maxamp(treeind(t),1)) < e
                        ready(treeind(t)) = true;
                        if all(ready)
                            flag = true;
                        end
                    else
                        amp(treeind(t)) = (maxamp(treeind(t),1)-minamp(treeind(t),1))/2 + minamp(treeind(t),1);
                    end
                    fprintf('Reached Spike of spike\n');
                else
                    if counter == 1 || sflag(treeind(t))
                        amp(treeind(t)) = amp(treeind(t)) - step;
                    else
                        ostep(treeind(t)) = ostep(treeind(t)) / 2;
                        amp(treeind(t)) = amp(treeind(t)) - ostep(treeind(t));
                    end
                    fprintf('Reached Voltage: Spike\n');
                end
                continue
                
            end
            if counter == 1
                Vrest(treeind(t)) = mean(thisv);
            end
            sflag(treeind(t)) = 0;
             fprintf('Reached Voltage: %2.1f, ',out.record{t}.cell.v{1}(end));
            if ischar(desv) && strcmpi(desv,'Spike')
                if isnan(maxamp(treeind(t))) || isnan(minamp(treeind(t)))
                    if isnan(minamp(treeind(t)))
                        minamp(treeind(t),:) = [amp(treeind(t)) thisv];
                    end
                    ostep(treeind(t)) = +step;
                    amp(treeind(t)) = amp(treeind(t)) + ostep(treeind(t));
                elseif maxamp(treeind(t),2) > -40 && abs(minamp(treeind(t),1)-maxamp(treeind(t))) < e
                    ready(treeind(t)) = true;
                    if all(ready)
                        flag = true;
                    end
                else
                    minamp(treeind(t),:) = [amp(treeind(t)) thisv];
                    amp(treeind(t)) = (maxamp(treeind(t))-minamp(treeind(t)))/2 + minamp(treeind(t));
                end
                
            else
               
                if abs(thisv - desv) <= e
                    ready(treeind(t)) = true;
                    if all(ready)
                        flag = true;
                    end
                else
                    if isnan(maxamp(treeind(t))) || isnan(minamp(treeind(t)))
                        if desv < thisv
                            maxamp(treeind(t),:) = [amp(treeind(t)) thisv];
                        else
                            minamp(treeind(t),:) = [amp(treeind(t)) thisv];
                        end
                        if isnan(maxamp(treeind(t))) || isnan(minamp(treeind(t)))
                            ostep(treeind(t)) = sign(desv-thisv)*step;
                            amp(treeind(t)) = amp(treeind(t)) + ostep(treeind(t));
                        else
                            amp(treeind(t)) = (maxamp(treeind(t),1)-minamp(treeind(t),1))/2 + minamp(treeind(t),1);
                            
                        end
                    else
                        if desv < thisv
                            maxamp(treeind(t),:) = [amp(treeind(t)) thisv];
                        else
                            minamp(treeind(t),:) = [amp(treeind(t)) thisv];
                            
                        end
                        amp(treeind(t)) = (maxamp(treeind(t))-minamp(treeind(t)))/2 + minamp(treeind(t));
                    end
                end
            end
        end
        if ischar(desv) && strcmpi(desv,'Spike')
            fprintf(' of spike\n');
        else
            fprintf(' of target voltage %g mV\n',desv);
        end
    end
    
    [~,~,ib] = intersect(find(ready),treeind);
    tree(ib) = [];
    neuron.mech(ib) = [];
    neuron.record(ib) = [];
    neuron.pp(ib) = [];
    neuron.APCount(ib) = [];
    treeind(ib) = [];
end
if counter > counterthresh
    warndlg('Not all target voltages could be reached! Check!')
end
