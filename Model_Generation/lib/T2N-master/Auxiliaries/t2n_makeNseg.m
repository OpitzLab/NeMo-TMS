function minterf = t2n_makeNseg(tree,minterf,par,mech)
% This function calculates at which position of each cell there will be a 
% segment node in NEURON. This is necessary for T2N to correctly translate
% TREES toolbox node information into NEURON section and segment information.
%
% INPUTS
% tree      single TREES toolbox tree morphology
% minterf	Nx3 mapping matrix created by neuron_template_tree (in morphology
%           folder)
% par       T2N neuron.params structure
% mech      T2N mech structure that belongs to the 'tree' morphology
%
% OUTPUT
% minterf	updated mapping matrix with 4th column showing the real segment node locations
%
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************

dodlambda = 0;
doeach = 0;

if ischar(par.nseg)
    %     pl = Pvec_tree(tree);  % path length of tree..
    len = len_tree(tree);
    idpar = idpar_tree(tree);
    
    if strcmpi(par.nseg,'dlambda')
        dodlambda = 1;
        
        D =  tree.D;
    elseif ~isempty(strfind(par.nseg,'ach')) % check for "each"
        each = double(cell2mat(textscan(par.nseg,'%*s %d'))); % get number
        doeach = 1;
    end
    
end
if ~isfield(par,'freq')
    par.freq = 300;
end
if ~isfield(par,'dlambda')
    par.dlambda = 0.1;
end
minterf(:,4) = 0;

for sec = 0:max(minterf(:,2))  %go through all sections
    secstart = find(minterf(:,2) == sec & minterf(:,3) == 0,1,'first');
    secend = find(minterf(:,2) == sec & minterf(:,3) == 1,1,'last');
    if dodlambda || doeach
        secnodestart = minterf(secstart,1);
        if secnodestart == 0  % this means the current section is the tiny section added for neuron... this should have nseg = 1
            secnodestart = 1;
            flag = true;
        else
            flag = false;
        end
        secnodestart2 = minterf(secstart+1,1);
        secnodeend = minterf(secend,1);
        
        if dodlambda
            if isfield(tree,'rnames') && isfield(mech,tree.rnames{tree.R(secnodestart2)}) && isfield(mech.(tree.rnames{tree.R(secnodestart2)}),'pas') && all(isfield(mech.(tree.rnames{tree.R(secnodestart2)}).pas,{'Ra','cm'}))
                Ra = mech.(tree.rnames{tree.R(secnodestart2)}).pas.Ra;
                cm = mech.(tree.rnames{tree.R(secnodestart2)}).pas.cm;
            elseif isfield(mech,'all') && isfield(mech.all,'pas') && all(isfield(mech.all.pas,{'Ra','cm'}))
                Ra = mech.all.pas.Ra;
                cm = mech.all.pas.cm;
            else
                %NEURON standard values for Ra and cm
                if isfield(tree,'rnames')
                    warning('Ra or cm of region %s in tree %s not specified',tree.rnames{tree.R(secnodestart)},tree.name)
                else
                    warning('Cannot find passive parameters for nseg calculation! If this is desired, you should define a fixed nseg value')
                end
                Ra = 35.4;
                cm = 1;
            end
        end
        
        if flag
            L = 0.0001;   % this is the length according to root_tree
        else
            L = sum(len(secnodestart2:secnodeend)); %length of section
        end
        if dodlambda
            lambda_f = sum(len(secnodestart2:secnodeend)./sqrt([D(secnodestart2);D(idpar(secnodestart2+1:secnodeend))] + D(secnodestart2:secnodeend)));
            
            lambda_f = lambda_f * sqrt(2) * 1e-5*sqrt(4*pi*par.freq*Ra*cm);
            
            if lambda_f == 0
                lambda_f = 1;
            else
                lambda_f = L/lambda_f;
            end
            nseg = floor((L/(par.dlambda*lambda_f)+0.9)/2)*2 + 1;
        else
            nseg = floor(L/each)+1;
        end
    else
        nseg = par.nseg;
    end
    if isfield(par,'accuracy')
        if par.accuracy == 2 || (par.accuracy == 1 && (~isempty(strfind(tree.rnames(tree.R(minterf(secend))),'axon')) || ~isempty(strfind(tree.rnames(tree.R(minterf(secend))),'soma'))) ) %triple nseg if accuracy is necessary also finds axonh
            nseg = 3 * nseg;
        end
    end
    pos = (2 * (1:nseg) - 1) / (2*nseg);    %calculate positions
    fac = (secend-secstart+1)/nseg;
    if fac > 1   % more tree nodes than segments
        [~,ind] = min(abs(repmat(minterf(secstart+1:secend,3),1,numel(pos)) - repmat(pos,secend-secstart,1)),[],2);  %find position to node which is closest to next segment location
        minterf(secstart+1:secend,4) = pos(ind);
    else     % more segments than tree nodes, so one has to add entries in minterf
        dupl = 0;
        for p = 1:numel(pos)
            ind = find(min(unique(abs(minterf(secstart+1:secend+dupl,3) - pos(p)))) == abs(minterf(secstart+1:secend+dupl,3) - pos(p)),1,'last'); % might be more cpu-consuming but attaches next pos in right order
            if minterf(secstart+1+ind-1,4) ~= 0
                minterf = minterf([1:secstart+1+ind-1,secstart+1+ind-1:end],:); %duplicate entry
                minterf(secstart+1+ind,4) = pos(p);  %overwrite duplicated entry with new value
                dupl = dupl +1;
            else
                minterf(secstart+1+ind-1,4) = pos(p);                % hier evt ausnahme für anfang und ende der section (=0)
            end
        end
    end
    minterf(secstart,4) = NaN;
end
end

