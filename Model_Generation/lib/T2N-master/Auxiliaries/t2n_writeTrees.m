function tree = t2n_writeTrees(tree,tname,savepath,server)
% This function transforms the tree file into hoc code and also saves a 
% interface file for correct node-section assignment
%
% INPUTS
% tree              tree cell array with morphologies (see documentation)
% tname             (optional) name after which all tree files will be 
%                   named (with counting). If 	not provided, the tree.name
%                   will be used, or, if not given simply 	“Tree”+number
% savepath          (optional) if this file destination string is given, 
%                   the function does not 	have to ask for the tree file 
%                   name to save via dialog
%
%OUTPUT
% tree              tree cell array with each tree containing a unique NEURON ID (NID)
% server            (optional for cluster mode, not fully implemented) server structure
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

options = '-w';


sflag = 0;
if isstruct(tree)
    tree = {tree};
    sflag = 1;
end
if size(tree,1) ~= numel(tree)  % check for correct 1st dimension 
    tree = tree';
end
morphfolder = fullfile(pwd,'morphos','hocs');

if exist('server','var')
    nrn_morphfolder = fullfile(server.clpath,'morphos/hocs');

    sshfrommatlabissue(server.connect,sprintf('mkdir -p %s',nrn_morphfolder));
else
    nrn_morphfolder = morphfolder;
end
nrn_morphfolder = regexprep(nrn_morphfolder,'\\','/');

if ~exist(morphfolder,'dir')
    mkdir(morphfolder);
end

orderchanged = false;
badchars = 0;

artflag = cellfun(@(x) isfield(x,'artificial'),tree);  % get boolean for trees being artificial
indartflag = find(artflag);   % get indices
[~,ia] = unique(cellfun(@(x) x.artificial,tree(artflag),'UniformOutput',0));  % find artificial cells that are of the same type
indWrite = cat(1,find(~artflag),indartflag(ia));  % only write trees that are not artificial and one artificial tree of each type
% tree = cat(1,tree(~artflag),tree(indartflag(ia)));  % only write trees that are not artificial and one artificial tree of each type

tim = tic;
if ~isempty(strfind(options,'-w'))
    w = waitbar(0,'Trees are transformed to hoc, please wait...');
end
for t=1:numel(tree)     % make neuron templates from trees and save/get minterface file
    if ~artflag(t)
        [tree{t}, order] = sort_tree(tree{t},'-LO');
        if ~all(order==sort(order))
            orderchanged = true;
        end
    end
    nonameflag = false;
    if ~artflag(t) && exist('tname','var') && ischar(tname) && ~isempty(tname)
        treename = tname;
        countupflag = true;
    elseif artflag(t) && ~isfield(tree{t},'name')
        treename = tree{t}.artificial;
        countupflag = false;
        nonameflag = true;
    elseif isfield(tree{t},'name') % this is also true if name exist and it is artificial
        treename = tree{t}.name;
        countupflag = false;
    else
        treename = 'Tree';
        countupflag = true;
        nonameflag = true;
    end
    if any(strfind(treename,'%'))
        badchars = badchars +numel(strfind(treename,'%'));
        treename(strfind(treename,'%')) = [];
    end
    if any(strfind(treename,'-'))
        badchars = badchars +numel(strfind(treename,'%'));
        treename(strfind(treename,'-')) = '_';
    end
    if any(strfind(treename,'.'))
        badchars = badchars +numel(strfind(treename,'.'));
        treename(strfind(treename,'.')) = '_';
    end
    if (numel(treename) < 5 || ~strcmp(treename(1:5),'cell_'))
        treename = strcat('cell_',treename);
    end
    if countupflag
        treename = sprintf('%s_%d',treename,t);
    end
%     if strfind(options,'-cl')
%         [server.connect, answer] = sshfrommatlabissue(server.connect,sprintf('ls %s/%s.hoc',nrn_morphfolder,treename));
%         fchk =  ~isempty(answer{1});
%     else
%         fchk = exist(fullfile(morphfolder,sprintf('%s.hoc',treename)),'file');
%     end
if any(t == indWrite)  % inly rewrite artificial trees once
    oname = treename;
    neuron_template_tree (tree{t}, fullfile(morphfolder,sprintf('%s.hoc',treename)), '-m');
    
    if exist('server','var')   %transfer files to server
        server.connect = sftpfrommatlab(server.connect,fullfile(morphfolder,sprintf('%s.hoc',oname)),sprintf('%s/%s.hoc',nrn_morphfolder,oname));
        pause(0.1)
        server.connect = sftpfrommatlab(server.connect,fullfile(morphfolder,sprintf('%s_minterf.dat',oname)),sprintf('%s/%s_minterf.dat',nrn_morphfolder,oname));
        pause(0.1)
        server.connect = sftpfrommatlab(server.connect,fullfile(morphfolder,sprintf('%s_minterf.mat',oname)),sprintf('%s/%s_minterf.mat',nrn_morphfolder,oname));
    end
end
    tree{t}.NID = treename;
    if nonameflag || badchars > 0
        tree{t}.name = treename;
    end
    if ~isempty(strfind(options,'-w'))
        waitbar(t/numel(tree),w);
    end
end

if sflag
    tree = tree{1};
end
if ~all(artflag)
    if nargin < 3
        display('Please resave trees to have the NEURON ID in each tree')
        save_tree(tree);
    else
        save_tree(tree,savepath);
    end
end
if ~isempty(strfind(options,'-w'))
    close(w)
end
if badchars > 0
    warning('Caution! %d bad chars had to be removed or replaced from the tree names since they cause writing errors! Please be sure to not use "%%" and "." in the names',badchars);
end
    tim = toc(tim);
    fprintf(sprintf('Tree hoc writing time: %g min %.2f sec\n',floor(tim/60),rem(tim,60)))

if orderchanged && nargout == 0
    warning('Caution, the node order of some trees had to be changed! Sort your trees with "sort_tree" to obtain the correct results')
end
