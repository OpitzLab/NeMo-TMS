function t2n_bAP(neuron,tree,amp,targetfolder_data,simplicity)
% This function performs an experiment on backpropagating action potentials 
% (bAPs) by "zapping" the cell(s) with a short high current and recording the
% membrane potential at all locations of the cell(s). The results are saved
% in a mat file named Exp_bap and the specification in neuron.experiment
% and can then be analyzed with t2n_plotbAP.
%
% INPUTS
% neuron            t2n neuron structure with already defined mechanisms
% tree              tree cell array with morphologies
% amp               amplitude [nA] of zap. Default is 1.3 nA
% targetfolder_data destination of temporary results file
% simplicity        0 (DEFAULT) measure voltage over whole dendritic tree
%                   1 only record from every second node
%                   2 only record voltages from the first primary dendrite
% 
% *****************************************************************************************************
% * This function is part of the T2N software package.                                                *
% * Copyright 2016-2019 Marcel Beining <marcel.beining@gmail.com>                                    *
% *****************************************************************************************************


neuron.params.v_init = -85.4;
if ~exist('amp','var')
    amp = 1.3; % nA
end
if ~exist('simplicity','var')
    simplicity = 0;
end

neuron.params.tstop = 1000;
neuron.params.dt=0.025;
neuron.params.cvode = 1;
nodes = cell(numel(tree),1);
plen = nodes;
eucl = nodes;
plotvals = nodes;  
bAP = nodes;
ipar = nodes;

hstep = t2n_findCurr(neuron,tree,neuron.params.v_init); %assuming a HP of xxx mV


for t = 1:numel(tree)
    
    plen{t} = Pvec_tree(tree{t});
    ipar{t} = ipar_tree(tree{t});
    ipar{t} = ipar{t}(T_tree(tree{t}),:);  % only paths from termination points
    ipar{t}(ipar{t}==0) = 1;
    
    if simplicity >= 2
        nodes{t} = unique(ipar(1,:));
    else
        nodes{t} = unique(ipar{t});
    end
    
    % this part would have been done by t2n anyway, however to avoid
    % loading a lot of redundant values into Matlab, nodes are reduced to
    % the locations were NEURON actually calculates voltage here
    minterf = load(fullfile(pwd,'morphos','hocs',sprintf('%s_minterf.mat',tree{t}.NID)));
    minterf = t2n_makeNseg(tree{t},minterf.minterf,neuron.params,neuron.mech{t});
    inode = zeros(numel(nodes{t}),1);
    for in = 1:numel(nodes{t})
        inode(in) = find(minterf(:,1) == nodes{t}(in),1,'first');    %find the index of the node in minterf
    end
    [~,ia] = unique(minterf(inode,[2,4]),'rows');
    nodes{t} = sort(nodes{t}(ia));
    if simplicity == 1   % reduce number of real recorded nodes to every second node.
        nodes{t} = nodes{t}(1:3:end);
    end
    
    neuron.record{t}.cell = struct('node',nodes{t},'record','v');
    neuron.pp{t}.IClamp = struct('node',1,'times',[-200 30,32.5],'amp', [hstep(t) hstep(t)+amp hstep(t)]); %n,del,dur,amp
    eucl{t} = eucl_tree(tree{t});
end
[out, ~] = t2n(neuron,tree,'-w-q-d');
tim = out.t;
for t = 1:numel(tree)
    plotvals{t} = NaN(numel(tree{t}.X),1);
    for x = 1:numel(nodes{t})
        [mx, ind] = max(out.record{t}.cell.v{nodes{t}(x)});
        basl = mean(out.record{t}.cell.v{nodes{t}(x)}(tim>=0 & tim<=30));
        ind2 = find(out.record{t}.cell.v{nodes{t}(x)} > (mx-basl)/2 + basl,1,'first');
        
        bAP{t}(x,:) = [nodes{t}(x) tim(ind) plen{t}(nodes{t}(x)) eucl{t}(nodes{t}(x)) mx basl tim(ind2)]; % ind nodes, time of max amp, PL at nodes, eucl at nodes, max amplitude, baseline, time of halfmax amp
        plotvals{t}(nodes{t}(x)) = mx;
    end
    % as not all nodes were recorded, interpolate the value for the nodes
    % between the recorded ones (only affects plotting the tree, not the
    % data graphs)
    for x = 1:size(ipar{t},1)
        plotvals{t}(ipar{t}(x,:)) = interp1(plen{t}(intersect(nodes{t},ipar{t}(x,:))),plotvals{t}(intersect(nodes{t},ipar{t}(x,:))),plen{t}(ipar{t}(x,:)),'pchip');
    end
    
end

save(fullfile(targetfolder_data,sprintf('Exp_bAP_%s.mat',neuron.experiment)),'bAP','plotvals','nodes','neuron','tree','tim')
