function [imp,freq] = t2n_resonance(neuron,tree,amp,holding_voltage,doerrorbar)
% Perform a resonance test (impedance measurement using oscillating current
% injections) on each cell and plot the result.
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms
% tree              tree cell array with morphologies (see documentation)
% amp               amplitude [nA] of the oscillating current injection
% holding_voltage   holding potential [mV] before current injection
% doerrorbar        (optional) Boolean if error bar should be added to plot
%
% OUTPUTS
% imp               matrix with impedances [MOhm] of each cell at all frequencies
% freq              frequency vector [Hz] same size as imp
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if ~exist(doerrorbar,'var')
    doerrorbar = 0;
end
if ~exist(holding_voltage,'var')
    holding_voltage = -70;
end

tvec=0:neuron.params.dt:neuron.params.tstop;
freq = (0.5*(tvec+neuron.params.dt)/(1000));
vec = sin(2*pi*tvec.*freq/1000) * amp;
figure;plot(tvec/1000,freq)
xlabel('time [s]')
ylabel('frequency [Hz]')


hold all;
xlabel('Frequency [Hz]')
ylabel('Impedance [M\Omega]')



hstep = t2n_findCurr(neuron,tree,holding_voltage,[],'-q-d');
for t = 1:numel(tree)
    neuron.pp{t}.IClamp = struct('node',1,'times',-200,'amp',hstep(t)); %n,del,dur,amp
    neuron.play{t}.IClamp = struct('node',1,'play','amp','times',tvec,'value',hstep(t)+vec); %n,del,dur,amp
    neuron.record{t}.cell = struct('node',1,'record','v');
end
% neuron_orig = neuron;
out = t2n(neuron,tree,'-q-d-w');
if any(cellfun(@(x) any(x.cell.v{1}>-30),out.record))
    warning('Caution! Spike was elicited!')
end


for t = 1:numel(tree)
    [pks,locs] = findpeaks(out.record{t}.cell.v{1},'MinPeakHeight',holding_voltage);
    [pks2,locs2] = findpeaks(-out.record{t}.cell.v{1},'MinPeakHeight',holding_voltage);
    [locs,ia] = sort([locs;locs2]);
    pks = [pks;-pks2];
    pks = pks(ia)-holding_voltage;
    imp(:,t) = abs(pks/amp);
end
if doerrorbar
    errorbar(freq(locs(2:2:end)-1),mean(imp(2:2:end,:),2),std(imp(2:2:end,:),[],2)/sqrt(numel(tree)),'color','k','linestyle','--');
    errorbar(freq(locs(1:2:end)-1),mean(imp(1:2:end,:),2),std(imp(1:2:end,:),[],2)/sqrt(numel(tree)),'color','k');
else
    plot(freq(locs(2:2:end)-1),mean(imp(2:2:end,:),2),'color','k','linestyle','--');
end
