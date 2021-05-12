function fig = t2n_plotVoltSteps(targetfolder_data,neuron,steps,subtract_hv)
% This function plots one or multiple voltage steps, which had been 
% previously simulated and saved with t2n_voltSteps. The exact loaded simulation
% is defined by neuron.experiment.
%
% INPUTS
% targetfolder_data     destination of temporary results file that was 
%                       provided to t2n_voltSteps
% neuron                t2n neuron structure with already defined mechanisms
% steps                 (optional) make another plot where only the
%                       specific current steps [nA] are plotted
% subtract_hv           Boolean if the current at the holding potential
%                       (before the step) should be subtracted. Default is 0
%
% OUTPUT
% fig                   figure handles to the plotted figures
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

loadingfile = t2n_catName(targetfolder_data,'Exp_VoltSteps',neuron.experiment,'.mat');
if ~exist('steps','var') || isempty(steps)
    steps = -120-12.1;
end

if ~exist('subtract_hv','var') || isempty(subtract_hv)
    subtract_hv = 0;
end

load(loadingfile,'dur','currVec','vstepsModel','tree')

fig(1) = figure; hold all
fig(2) = figure;hold all
if subtract_hv
    for f = 1:numel(tree)
        for s = 1:size(currVec,2)
            currVec{f,s}(2,:) = currVec{f,s}(2,:) - mean(currVec{f,s}(2,currVec{f,s}(1,:) >=0 & currVec{f,s}(1,:) < dur(1))); 
        end
    end
end

str =  '';
figure(fig(1))
for f = 1:numel(tree)
    for s = 1:size(currVec,2)
        subplot(floor(sqrt(size(currVec,2))),ceil(sqrt(size(currVec,2))),s)
        hold all
        plot(currVec{f,s}(1,:),currVec{f,s}(2,:),'LineWidth',3,'LineStyle','-','Color',tree{f}.col{1})
        ylabel('Current [pA]')
        xlabel('Time [ms]')
        ylim([-200,200])
        title(sprintf('VClamp % 4.4g mV%s',vstepsModel(s),str));
        if any(vstepsModel(s) == steps)
            figure(fig(2))
            plot(currVec{f,s}(1,:),currVec{f,s}(2,:),'LineWidth',1,'LineStyle','-','Color',tree{f}.col{1})
            ylabel('Current [pA]')
            xlabel('Time [ms]')
            ylim([-500,200])
            title(sprintf('VClamp % 4.4g mV%s',vstepsModel(s),str));
            xlim([0 300])
            figure(fig(1))
        end
    end
end