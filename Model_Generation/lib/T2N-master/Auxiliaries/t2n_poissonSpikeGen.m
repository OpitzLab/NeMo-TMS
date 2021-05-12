function [spikeMat, tVec] = t2n_poissonSpikeGen(freq, par, nTrials)
% This function creates Poisson spike trains at a frequency 'freq' the
% output is a spike matrix and a time vector, both which can be used to
% define a artificial NEURON VecStim that follows such a spike train.
%
% INPUTS
% freq          desired frequency of the spiking [Hz]
% par           parameter structure of t2n neuron.params (see documentation)
% nTrials       number of independent spike traces to be generated ( =rows
%               in the spike matrix)
%
% OUTPUTS
% spikeMat      logical matrix with ones where a spike occurs, e.g. for
%               input to a NEURON VecStim point process
% tVec          corresponding time vector [ms]
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


if nargin < 3
    nTrials = 1;
end

nBins = floor(par.tstop/par.dt);
spikeMat = rand(nTrials, nBins) < freq*par.dt/1000;
tVec = 0:par.dt:par.tstop-par.dt;