function [currVec,out] = t2n_voltSteps(neuron,tree,vstepsModel,dur,holding_voltage,targetfolder_data)
% This function performs one or multiple voltage steps (in the squared
% pulse format, i.e. voltage step is surrounded by step to the baseline
% potential) in the cells given by "tree" and "neuron" and saves the 
% results in a mat file named according to neuron.experiment for later analysis.
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms
% tree              tree cell array with morphologies
% vstepsModel       vector with voltages [mV] that will be held
% dur               duration time period of voltage step [ms]. Can be a 1x3
%                   vector which then includes the pre and post step 
%                   duration at which the cell is held at the holding 
%                   potential, or a scalar with only the step duration and 
%                   pre/post being set to 100 ms
% holding_voltage   potential [mV] at which cell is held before and after 
%                   the voltage step
% targetfolder_data destination of temporary results file
%
% OUTPUTS
% currVec           cell array of vectors containing the simulation time 
%                   vector and the measured current during the simulation 
%                   for each cell (first dim) and voltage step (second dim)
% out               the direct output structure of t2n, necessary for some
%                   other t2n functions
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


if ~exist('holding_voltage','var')
    holding_voltage = -80;
end
if numel(holding_voltage) < 2
    holding_voltage = repmat(holding_voltage,1,2);
end
if ~exist('dur','var')
    dur = [100 100 100];
elseif numel(dur) == 1
    dur = [100 dur 100];
end

elecnode = 1;

neuron.params.tstop = sum(dur);


nneuron = cell(numel(vstepsModel),1);
nneuron{1} = neuron;

for s = 1:numel(vstepsModel)
    amp = cat(2,holding_voltage(1), vstepsModel(s), holding_voltage(2));
    for t = 1:numel(tree)
        nneuron{s}.pp{t}.SEClamp = struct('node',elecnode,'rs',15,'dur', dur,'amp', amp);
        if s == 1
            nneuron{s}.record{t}.SEClamp = struct('record','i','node',elecnode);
        end
    end
end
nneuron = t2n_as(1,nneuron);
out = t2n(nneuron,tree,'-q-w');
if any(cellfun(@(x) x.error,out(cellfun(@(x) isfield(x,'error'),out))))
    return
end

ind = vstepsModel == holding_voltage(1);
if ~any(ind)
   mholding_current = NaN(numel(tree),1); 
else
    for t = 1:numel(tree)
        mholding_current(t) = mean(out{ind}.record{t}.SEClamp.i{1}(find(out{ind}.t>=dur(1)+dur(2)*0.9,1,'first'):find(out{ind}.t<sum(dur(1:2)),1,'last')) *1000);  % get steady state voltage (electrode current at 90-100% of step duration)
    end
end
for s = 1:numel(vstepsModel)
    for t = 1:numel(tree)
        steadyStateCurrVec(s,t) =  mean(out{s}.record{t}.SEClamp.i{1}(find(out{s}.t>=dur(1)+dur(2)*0.9,1,'first'):find(out{s}.t<sum(dur(1:2)),1,'last')) *1000); % get steady state voltage (electrode current at 90-100% of step duration)
        currVec{t,s} =  [out{s}.t';out{s}.record{t}.SEClamp.i{1}' *1000];
    end
end
if nargout == 0
    save(fullfile(targetfolder_data,sprintf('Exp_VoltSteps_%s.mat',neuron.experiment)),'dur','mholding_current','neuron','holding_voltage','steadyStateCurrVec','currVec','vstepsModel','tree')
end