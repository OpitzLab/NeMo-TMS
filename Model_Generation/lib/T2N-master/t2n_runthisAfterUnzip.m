% Run this code when you extracted the T2N zip file including this script at its final loation.
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

p = fileparts(mfilename('fullpath')); % get path to this script
addpath(genpath(p)); % add folder including subfolders to the Matlab path
fprintf('Added %s and its subdirectories to the Matlab path!\n',p);

if ispc
    lookForCommand = 'where';
    sep = '&&';
else
    lookForCommand = 'which';
    sep = ';';
end
[~,outp] = system([lookForCommand, ' git']);
% check if any git has been found
if isempty(outp) || ~isempty(regexp(outp,'not found','ONCE'))  || ~isempty(regexp(outp,'no nrniv','ONCE')) || ~isempty(regexp(outp,'not find','ONCE')) || strcmp(outp,sprintf('\n')) || isempty(regexp(outp,'git','ONCE'))
    error('Git was not found on your system, thus masterplotter submodule could not be installed. Please install git (e.g. https://gitforwindows.org/, activate option adding git to PATH) and rerun this script file!')
else
    [~,outp] = system(sprintf('cd "%s" %s git submodule update --init --recursive',p,sep));
end
