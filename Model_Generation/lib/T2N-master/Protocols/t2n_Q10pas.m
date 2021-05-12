function neuron = t2n_Q10pas(neuron,celsius)
% This function adjusts the passive parameters g and Ra to a given temperature
% (assuming that g and Ra were defined for 24°C before). Hence do not use it
% multiple times!
%
% INPUTS
% neuron                t2n neuron structure with already defined mechanisms
% celsius               the temperature [C] to which the passive parameters
%                       will be adjusted to
% 
% OUTPUT
% neuron                the modified neuron structure
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

scaleg = 1.98^((celsius-24)/10);
scaleRa = 0.8^((celsius-24)/10);

for t = 1:numel(neuron.mech)
    fields = fieldnames(neuron.mech{t});
    for f1 = 1:numel(fields)
        if isfield(neuron.mech{t}.(fields{f1}),'pas')
            neuron.mech{t}.(fields{f1}).pas.g = neuron.mech{t}.(fields{f1}).pas.g * scaleg;
            neuron.mech{t}.(fields{f1}).pas.Ra = neuron.mech{t}.(fields{f1}).pas.Ra * scaleRa;
        end
    end
end

if isfield(neuron,'experiment') && ~isempty(neuron.experiment)
    neuron.experiment = strcat(neuron.experiment,'_Q10pas');
else
    neuron.experiment = 'Q10pas';
end
