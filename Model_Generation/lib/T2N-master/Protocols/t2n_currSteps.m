function t2n_currSteps(neuron,tree,targetfolder_data,ostruct)
% This function performs one or multiple current steps in the cells given
% by "tree" and "neuron" and saves the results in a mat file named
% according to neuron.experiment.
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms
% tree              tree cell array with morphologies
% targetfolder_data destination of temporary results file
% ostruct           structure with fields defining the curr step simulation
%                   amp     vector with amplitudes [nA]
%                   delay   time point at which current injection starts [ms]
%                   duration time period of current injection [ms]
%                   holding_voltage (optional) potential at which cell is held before current injection [mV]
%                   spikeThresh (optional) threshold above which spikes are
%                                          detected (default is -15 mV)
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 4
    ostruct = struct();
end
if isfield(ostruct,'recordnode') && isnumeric(ostruct.recordnode) && numel(ostruct.recordnode) == numel(tree)
    recordnode = ostruct.recordnode;
else
    recordnode = ones(numel(tree),1);
end
if isfield(ostruct,'stimnode') && isnumeric(ostruct.stimnode) && numel(ostruct.stimnode) == numel(tree)
    stimnode = ostruct.stimnode;
else
    stimnode = ones(numel(tree),1);
end
if ~isfield(ostruct,'duration')
    ostruct.duration = 200;  % standard duration 200 ms
end
if ~isfield(ostruct,'numAP')
    ostruct.numAP = 0;
end
if ~isfield(ostruct,'spikeThresh')
    ostruct.spikeThresh = -10;  % std spike thresh of -10 mV
end
if ~isfield(ostruct,'amp')
    ostruct.amp = (0:5:90)/1000;  % standard current steps 0-90 pA
end
if ~isfield(ostruct,'delay')
    ostruct.delay = 55.5;  % standard current steps 0-90 pA
end

neuron.params.accuracy = 1;  % for more nseg in axon and soma!

if isfield(ostruct,'coarse') && ostruct.coarse == 1
    neuron.params.nseg = 1;
    neuron.params.dt=0.1;  % does only count if cvode = 0
elseif isfield(ostruct,'coarse') && ostruct.coarse == 0.5
    neuron.params.dt=0.05;  % does only count if cvode = 0
else
    neuron.params.dt=0.025;  % does only count if cvode = 0
end
cstepsSpikingModel = ostruct.amp;  % 0:5:120

neuron.params.tstop = 150+ostruct.delay+ostruct.duration;

if isfield(ostruct,'holding_voltage') && ~isnan(ostruct.holding_voltage)
    hstep = t2n_findCurr(neuron,tree,ostruct.holding_voltage,[],'-q-d');
else
    hstep = zeros(1,numel(tree));
end

if isnan(hstep)
    return
end
for t=1:numel(tree)
    neuron.APCount{t} = [1,ostruct.spikeThresh];
end

nneuron = cell(numel(cstepsSpikingModel),1);
nneuron{1} = neuron;
for s = 1:numel(cstepsSpikingModel)
    for t = 1:numel(tree)
        nneuron{s}.pp{t}.IClamp = struct('node',stimnode(t),'times',[-100,ostruct.delay,ostruct.delay+ostruct.duration],'amp', [hstep(t) hstep(t)+cstepsSpikingModel(s) hstep(t)]); %n,del,dur,amp
        nneuron{s}.record{t}.cell = struct('node',recordnode(t),'record','v');
    end    
end
nneuron = t2n_as(1,nneuron);

if ostruct.numAP > 0
    amp = t2n_findFreq(nneuron{1},tree,ostruct.numAP,'-q-d');
    for t = 1:numel(tree)
        nneuron{1}.pp{t}.IClamp.amp = [hstep(t) amp(t) hstep(t)]; %n,del,dur,amp  %WICHTIG! nur amp da hstep nicht abgezogen
    end
end

out = t2n(nneuron,tree,'-q-d-w'); % run simulations

numspikes = zeros(numel(tree),numel(cstepsSpikingModel));
voltVec = cell(numel(tree),numel(cstepsSpikingModel));
timeVec = voltVec;

for s = 1:numel(cstepsSpikingModel)
    for t = 1:numel(tree)
        if isfield(out{s},'error') && out{s}.error > 0
            voltVec{t,s} = [] ;
            timeVec{t,s} = [];
            numspikes(t,s) = NaN;
        else
            voltVec{t,s} = out{s}.record{t}.cell.v{1} ;
            timeVec{t,s} = out{s}.t;
            numspikes(t,s) = numel(out{s}.APCtimes{t}{1});
        end
    end
end
simDef = ostruct;
save(fullfile(targetfolder_data,sprintf('Exp_Spiking_%s.mat',neuron.experiment)),'voltVec','timeVec','numspikes','cstepsSpikingModel','tree','nneuron','simDef')
