function t2n_plotRaster(varargin)
% This function plots the raster plot of one or multiple spike trains that 
% were obtained from an APCount, created manually or generated e.g. with 
% t2n_poissonSpikeGen. 
%
% INPUTS
% spikeMat      logical matrix with Ones where a spike occurs
%               alternatively it can be a vector with spiking times, or for
%               multiple spike trains, a cell array with each cell comprising
%               the spike time vectors
% tVec          Corresponding time vector [ms]. Only necessary if spikeMat
%               is a logical matrix
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


switch nargin
    case 1
        spikeMat = varargin{1};
    case 2
        spikeMat = varargin{1};
        tVec = varargin{2};
end
if ~islogical(spikeMat) % if spikeMat is not spikeMat but times of spike
    if ~exist('tvec','var')
        if iscell(spikeMat)
            spikeMat = cellfun(@(x) x(:),spikeMat,'uni',0); % sort spike times in one direction
            tVec = sort(cat(1,spikeMat{:})); % use all sorted spike times as tvec
        else
            tVec = sort(spikeMat);  % use all sorted spike times as tvec
        end
    else
        error('time vector has to be provided if spikeMat is a logical matrix')
    end
        
    tmp = spikeMat;
    if ~iscell(tmp)
        tmp = {tmp};
    end
    spikeMat = false(numel(tmp),numel(tVec)); % init spikeMat
    for n = 1:numel(tmp)
        if iscell(tmp{n})
            tmp{n} = tmp{n}{1};
        end
        spikeMat(n,interp1(tVec,1:numel(tVec),tmp{n})) = true; %make spikeMat 1 at times
    end
end

hold all;
for trialCount = 1:size(spikeMat,1)
    spikePos = tVec(spikeMat(trialCount, :));
    for spikeCount = 1:length(spikePos)
        plot([spikePos(spikeCount) spikePos(spikeCount)], ...
            [trialCount-0.4 trialCount+0.4], 'k');
    end
end
ylim([0 size(spikeMat, 1)+1]);
xlim([0 tVec(end)])
set(gca,'YTick',(1:trialCount))
set(gca,'YTickLabel',1:1:trialCount)
ylabel('Cell number')
xlabel('Time [ms]')