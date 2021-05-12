function t2n_renameNrnmech(newname,path)
% This function renames the nrnmech dll file to the name specified in newname. 
% Also deletes all .o and .c files created during dll compilation by
% NEURON.
%
% INPUTS
% newname	string with new name for dll file
% path      (optional) path to the main folder of the model (containing the
%           folder lib_mech). If not provided, the current Matlab working
%           directory is used.
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


if nargin < 1
    newname = 'nrnmech.dll';
end
if nargin < 2
    path = pwd;
end
origpth = pwd;

if isnumeric(newname)
    nam = sprintf('nrnmech_win%d.dll',newname);
elseif strcmp(newname(end-3:end),'.dll')
    nam = newname;
else
    nam = sprintf('%s.dll',newname);
end

if isempty(strfind(path,'lib_mech'))
    if ~exist(fullfile(path,'lib_mech'),'file')
        error('No folder lib_mech exists')
    else
       cd(fullfile(path,'lib_mech'))
    end
end
if ~strcmp('nrnmech.dll',nam)
    delete(nam)
    movefile('nrnmech.dll',nam)
end
fils = dir();
fils = fils(3:end);
fils = fils(cellfun(@(x) strcmp(x(end-1:end),'.o')|strcmp(x(end-1:end),'.c'),{fils.name}));
delete(fils.name)
cd(origpth)
display(sprintf('nrnmech.dll successfully renamed in %s ... and .o/.c files deleted',nam))