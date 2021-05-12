function n = t2n_getref(n,neuron,field)
% This function returns the neuron instance that is referenced in neuron
% instance 'n' at field 'field'. If the neuron instance has own definitions 
% at 'field', n itsself is returned. The function t2n_as produces such references.
%
% INPUTS
% n              	index of the neuron instance, which reference is
%                   searched for
% neuron            t2n neuron structure with already defined mechanisms (see documentation)
% field           	the field for which reference is searched for, e.g. params, mech, pp etc.
%
% OUTPUTS
% n                 reference to the neuron instance comprising the
%                   definitions of 'field'
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


if isfield(neuron{n},field)
    if strcmp(field,'tree')
        if isnumeric(neuron{n}.tree)        % n is already correct ref
            return
        elseif isempty(strfind(neuron{n}.tree,'sim'))
            n = NaN;
            return
        end
        
        while 1
            if ~isnumeric(neuron{n}.tree)  % ref to a sim
                if str2double(neuron{n}.tree(4:end)) == n   % ref to itsself
                    if n == 1
                        n = [];
                        return
                    else
                        n = NaN;
                        break
                    end
                else                        % ref to another sim
                    n = str2double(neuron{n}.tree(4:end));
                end
            else                        % it found the ref, break
                break
            end
        end
        
    else
        while 1
            if isfield(neuron{n},field)
                if isnumeric(neuron{n}.(field))  % ref to a sim
                    if neuron{n}.(field) == n   % ref to itsself
                        n = NaN;
                        break
                    else                        % ref to another sim
                        n = neuron{n}.(field);
                    end
                else                        % it found the ref, break
                    break
                end
            else  % ref to a sim which not that field
                n = NaN;
                break
            end
        end
    end
else
    n = NaN;
end

end