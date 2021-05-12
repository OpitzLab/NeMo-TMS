function t2n_plotCurrSteps(targetfolder_data,neuron,steps)
% This function plots one or multiple current steps, which had been 
% previously simulated and saved with t2n_currSteps. The exact loaded simulation
% is defined by neuron.experiment.
%
% INPUTS
% targetfolder_data     destination of temporary results file
% neuron                t2n neuron structure with already defined mechanisms
% steps                 (optional) restrict plot on specific current steps [nA]
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

load(t2n_catName(targetfolder_data,'Exp_Spiking',neuron.experiment,'.mat'))

if nargin < 3 || isempty(steps)
    steps = [];
end
hold all 
if any(cstepsSpikingModel==0.01)
    Rin = (cellfun(@(x) max(x),voltVec(:,cstepsSpikingModel==0.01))-cellfun(@(x) x(1),voltVec(:,cstepsSpikingModel==0.01)))/0.01;
    fprintf('Mean Rin in Model(@+10pA) is %g +- %g MOhm (s.e.m., -10mV)\n',mean(Rin),std(Rin)/sqrt(numel(Rin)))
end
for f=1:size(voltVec,1)
    for s = 1:size(voltVec,2)
        if isempty(steps)
            subplot(round(sqrt(size(voltVec,2))),ceil(sqrt(size(voltVec,2))),s)
            hold all
            if ~isempty(timeVec{f,s})
                plot(timeVec{f,s},squeeze(voltVec{f,s}),'LineWidth',1,'Color',tree{f}.col{1})
            end
        elseif any(cstepsSpikingModel(s) == steps)
            if ~isempty(timeVec{f,s})
                plot(timeVec{f,s},squeeze(voltVec{f,s}),'LineWidth',1,'Color',tree{f}.col{1})
            end
        end
        xlabel('Time [ms]')
        ylabel('Membrane voltage [mV]')
    end
end
