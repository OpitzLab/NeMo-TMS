function strct = t2n_catStruct(varargin)
% This function concatenates two structures which have the same fields.
%
% INPUTS
% varargin          several structures with same field names
%
% OUTPUT
% strct             concatenated structure
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

strct = struct();
for v = 1:nargin
    fields = fieldnames(varargin{v});
    for f = 1:numel(fields)
        if isfield(strct,fields{f})
            if ~isstruct(strct.(fields{f}))
                if (isnumeric(strct.(fields{f})) && ~all(strct.(fields{f}) == varargin{v}.(fields{f}))) || (ischar(strct.(fields{f})) && ~strcmp(strct.(fields{f}),varargin{v}.(fields{f}))) % check if values are not the same
                    warning('Entry of field %s (%d) has been overwritten with %d.\n',fields{f},strct.(fields{f}),varargin{v}.(fields{f}))
                end
                strct.(fields{f}) = varargin{v}.(fields{f});  % overwrite value
            else
                strct.(fields{f}) = t2n_catStruct(strct.(fields{f}),varargin{v}.(fields{f}));
            end
        else
            strct.(fields{f}) = varargin{v}.(fields{f});
        end
    end
end