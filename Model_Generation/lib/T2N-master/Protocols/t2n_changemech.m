function strct = t2n_changemech(strct,change,mode)
% This function changes parameters in the T2N neuron mech structure 
% according to argument "change".
%
% INPUTS
% strct     neuron structure containing ion channel densities or 
%           alternatively only the mech field of the neuron structure
% change    structure with field names according to mechanism parameter 
%           names that should be changed and values either describing the 
%           factor by which it should be changed or an absolute value
% mode      how parameters should be changed: 
%           1 or 'relative': (DEFAULT) values in change are relative factors (e.g. 0.5 for 50% decrease)
%           2 or 'absolute': values are absolute values that overwrite the values in strct
%
% OUTPUTS
% strct     updated neuron structure
%
% NOTE
% Cell arrays of neuron structure are not supported.
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

if nargin < 3 || isempty(mode) || (ischar(mode) && strcmpi(mode,'relative'))
    mode = 1;
else
    mode = 2;
end
changer = ones(1,2);  % this variable is used to either calculate a relative or the absolute value
flag = 0; % flag for neuron structure
if iscell(strct)
    error('Cell arrays of neuron structures are not supported as input');
end
if isfield(strct,'mech') % if whole neuron structure was the input, extract the mechanism structure and save the rest for later
    ostrct = strct;
    strct = strct.mech;
    flag = 1;
end
if ~isempty(change) && isstruct(change)
    % convert neuron parameter into mechanism name and parameter name
    chfield = fieldnames(change);
    chfieldCell = cellfun(@(x) strsplit(x,'_'),chfield,'uni',0);
    chfieldCell = cat(1,chfieldCell{:});
    % go through all neurons
    for t = 1:numel(strct)
        field = fieldnames(strct{t});  
        for f = 1:numel(field)  % go through all regions
            if isstruct(strct{t}.(field{f}))
                ffield = fieldnames(strct{t}.(field{f}));
                [~,~,ind] = intersect(ffield,chfieldCell(:,2));  % look for intersection between parameters to change and parameters in region
                if ~isempty(ind)
                    for in = 1:numel(ind)
                        changer(1) =  strct{t}.(field{f}).(chfieldCell{ind(in),2}).(chfieldCell{ind(in),1});
                        strct{t}.(field{f}).(chfieldCell{ind(in),2}).(chfieldCell{ind(in),1}) = changer(mode) * change.(chfield{ind(in)});
                    end
                end
%                 for ff = 1:numel(ffield)
%                     if isstruct(strct{t}.(field{f}).(ffield{ff}))
%                         ind = indect(fieldnames(strct{t}.(field{f}).(ffield{ff})),chfield);
%                         if ~isempty(ind)
%                             for in = 1:numel(ind)
%                                 changer(1) = strct{t}.(field{f}).(ffield{ff}).(chfield{ind(in)});
%                                 strct{t}.(field{f}).(ffield{ff}).(chfield{ind(in)}) = changer(mode) * change.(chfield{ind(in)});
%                             end
%                         end
%                     end
%                 end
            end
        end
    end
end
if flag  % put changed mech back into the neuron structure
    ostrct.mech = strct;
    strct = ostrct;
end