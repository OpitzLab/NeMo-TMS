function outname = t2n_catName(varargin)
% This function concatenates strings to a full file, e.g. for file names.
% If the first input is a path (i.e. containing separators / or \) than the
% function concatenates the strings to a full file path
% 
% INPUTS
% varargin      any number of strings
%
% OUTPUTS
% outname      concatenated full name
%
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************
Path = '';
if any(~isempty(strfind(varargin{1},'/'))) || any(~isempty(strfind(varargin{1},'\')))
    if any(~isempty(strfind(varargin{1},'/')))
        filesepa = '/';
    else
        filesepa = '\';
    end
    Path = varargin{1};
    varargin = varargin(2:end);
    if ~strcmp(Path(end),filesepa)
        Path(end+1) = filesepa;
    end
end
outname = '';
for n = 1:numel(varargin)
    if n == numel(varargin) && strcmp(varargin{n}(1),'.')
        outname = strcat(outname,varargin{n});
    else
        outname = strcat(outname,'_',varargin{n});
    end
end
outname = strcat(Path,outname(2:end));
