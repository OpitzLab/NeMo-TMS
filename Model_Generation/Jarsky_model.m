%% Jarsky CA1 PC model generator for TMS simulation
% Implements the Jarsky model of the CA1 PC:
% Jarsky T, Roxin A, Kath WL, Spruston N (2005) Conditional dendritic spike propagation following distal synaptic activation of hippocampal CA1 pyramidal neurons. Nat Neurosci 8:1667-76
% 
% Input: SWC or MTR file containing morphology data in ./morphos/
%  
% Adapted by Nicholas Hananeia, 2019-2021
%%
function[] = Jarsky_model()
fclose('all');


% install T2N and TREES and add path for necessary files
run('./lib/init_t2n_trees.m');

if ~exist('./Jarsky_files/lib_mech/nrnmech.dll', 'file') && ~exist('./Jarsky_files/lib_mech/x86_64/.libs/libnrnmech.so', 'file')
    msg = 'Compiled mod file not detected! Compile mods in Jarsky_files/lib_mech before continuing! Exiting...';
    msgbox(msg, 'Mods not compiled');
    return
end
    

addpath('Jarsky_files');
[input_cell, cell_path] = uigetfile({'*.swc';'*.mtr'}, 'Select morphology file', './morphos');

tstop                    = 1100;%40000;
dt                       = 0.05;

% Define standard parameters
neuron.params            = [];
neuron.params.celsius    = 35;
neuron.params.v_init     = -85;
neuron.params.prerun     = 200;
neuron.params.tstop      = tstop;
neuron.params.dt         = dt;
neuron.params.nseg       = 'dlambda';
neuron.params.dlambda    = 0.025;
neuron.params.freq       = 500;


fileid = strcat(cell_path, input_cell);
[~, name, ~] = fileparts(fileid);
trees = load_tree(fileid); %Specify input morphology here!
syn_distance = 100; %distance from soma at which synapse will be placed
resampling_interval = 0.5;

axon_type = menu('Choose desired axon:','Do not alter','No axon','Stick axon','Myelinated axon');

%%

for cell_num = 1:length(trees)
simName = inputdlg('Enter model name; leave blank to use cell filename', 'Model name');
%% Load morphology
treeFilename = './morphos/place_tree.mtr'; %Input file here!
treepath = '';

if length(simName{1})~=0 
    name = simName{1};
end

if(exist(strcat('../Models/', name), 'file'))
    msg = 'Model by that name already exists. Delete the model folder or use a different name. Exiting...';
    msgbox(msg, 'Duplicate model name');
    rmpath('Jarsky_files');
    return
end

opts.Interpreter = 'tex';
dist_input = inputdlg('Enter distance of synapse from soma on apical dendrite (\mum):', 'Enter synapse distance',1,{num2str(syn_distance)},opts);

if length(dist_input{1})~=0
    if str2double(dist_input{1}) >= 0  %we don't want a negative value causing nonsense
        syn_distance = str2double(dist_input{1});
    else
        msg = 'Negative distance entered. Reverting to default.';
        msgbox(msg, 'Warning: Negative distance entered.');
    end 
end

%to de-group the morphologies (if necessary), and have the different tree
%strucutres in the cell array 'tree':
current_dir = strcat('../Models/', name);
t2n_initModelfolders(current_dir);



%Load morphologies. 
if(length(trees) == 1)
    tree{1,1} = trees;
    neuron.params.exchfolder = strcat('../Models/',name, '/Code/');
else
    tree{1,1} = trees{cell_num};
    name = trees{cell_num}.name;
    neuron.params.exchfolder = strcat('../Models/',name, '/Code/');
end

if(axon_type == 2 || axon_type == 3)
    tree{1,1} = strip_axon(tree{1,1});
end


%% Divide the tree morphology (if it hasn't been divided before)
for t                    = 1 : numel (tree)
    % option j jarsky?, axon include axon
    %tree{t}.R            = tree{t}.R * 0 + 2;
    %tree{t}.R (1)        = 1;
    tree{t}.rnames       = {'soma', 'axon', 'dendrite' 'dendrite'};
    
    if(axon_type == 3)
        tree{1,1} = add_axon(tree{1,1});
    end
    if(axon_type == 4)
        tree{t} = myelinate_axon(tree{t});
    end
    tree{t}              = CA1pyramidalcell_sort_Jmodel_len(tree{t},'-j -axon');
end 


%% Convert the tree into NEURON
for t                    = 1 : numel (tree)
    if ~all (cellfun (@(x) isfield (x, 'NID'), tree)) || ...
            ~all (cellfun (@(x) exist (fullfile ( ...
            current_dir, 'morphos', 'hocs', [x.NID, '.hoc']), 'file'), ...
            tree))
        answer = 'OK';
        if strcmp        (answer, 'OK')
            tree{t}      = sort_tree      (tree{t}, '-LO');
            tname = strcat('Jarsky_',name);
            % Tanslation of morphologies into hoc file:
            tree         = t2n_writeTrees (tree,tname, fullfile (treepath, treeFilename));
            figure(cell_num);
            xplore_tree(tree{t}, '-2');
            xlabel('x(\mum)');
            ylabel('y(\mum)');
            zlabel('z(\mum)');
            title(name,'Interpreter','none');
        end
    end
end


%% Add passive parameters
cm                       = 0.75;              % Membrane capacitance (µF/cm²)
Ra                       = 200;               % Cytoplasmic resistivity (ohm*cm)
Rm                       = 40000;             % Membrane resistance (ohm/cm²) (uniform)
gpas                     = 1;
e_pas                    = -66;
for t                    = 1 : numel (tree)
    % do not scale spines:
    neuron.mech{t}.all.pas      = struct ( ...
        'cm',        cm,  ...
        'Ra',        Ra,  ...
        'g',         gpas / Rm,  ...
        'e',         e_pas);
end

%% Add active mechanisms
% To get the regions that should be ranged:
% taken the regions from tree 1 because all of them have the same region
% definition:
treeregions              = tree{1}.rnames;
noregions                = {'soma', 'basal', 'hill', 'iseg', 'myelin', 'node'};
x                        = false (size (treeregions));
for r                    = 1 : numel (noregions)
    % <-- Flag the ones that noregions{r} matches:
    x                    = x | ~cellfun (@isempty, strfind (treeregions, noregions{r}));
end
treeregions (x)          = [];    % <-- Delete all the flagged lines at once
% ********** Na conductance (gNabar)
nainfo.gbar              = 0.040;                % in S/cm2
nainfo.gnode             = 30.0;                % in S/cm2
nainfo.region            = treeregions;
% ********** Delayed rectifier K+ conducatnce (gKdr)
gkdr                     = 0.040;                       % in S/cm2 (uniform)
% ********** A-type K+ channel proximal (gAKp) and distal (gAKd)
kainfo.gka               = 0.048;                 % in S/cm2
kainfo.gka_ax            = kainfo.gka * 0.2;     % in S/cm2
kainfo.ek                = -77;
kainfo.region            = treeregions;
for t                    = 1 : numel (tree)
    % Distribution of the channels that depend on path distance
    vec_gNa{t}                    = range_conductanceNa (nainfo, tree{t}, ...
        '-wE'); % some option determining excitability
    vec_gKa{t}                    = range_conductanceKa (kainfo, tree{t});
    
    neuron.mech{t}.range.nax      = struct ( ...
        'gbar',                vec_gNa{t});
    neuron.mech{t}.all.nax        = struct ( ...
        'gbar',                nainfo.gbar, ...
        'ena',                 55);
    
    neuron.mech{t}.all.kdr        = struct ( ...
        'gkdrbar', gkdr, ...
        'ek',                  -77);
    
    neuron.mech{t}.kap            = struct ( ...
        'gkabar',              vec_gKa{t}.proximal);
    neuron.mech{t}.all.kap        = struct ( ...
        'gkabar',              kainfo.gka, ...
        'ek',                  -77);
    neuron.mech{t}.all.kad = struct();
    
    neuron.mech{t}.range.kad      = struct ( ...
        'gkabar',              vec_gKa{t}.distal);
    neuron.mech{t}.proxAp.kad     = struct ( ...
        'gkabar',              kainfo.gka, ...
        'ek',                  -77);
    neuron.mech{t}.middleAp.kad   = struct ( ...
        'gkabar',              kainfo.gka, ...
        'ek',                  -77);
    neuron.mech{t}.distalAp.kad   = struct ( ...
        'gkabar',              kainfo.gka, ...
        'ek',                  -77);
    neuron.mech{t}.tuft.kad       = struct ( ...
        'gkabar',              kainfo.gka, ...
        'ek',                  -77);
    %Myelin segments have lowered membrane capacitance
    neuron.mech{t}.myelin.pas = struct(...
        'cm', 0.01, ...
        'g_pas', 1/1.125e6);
    %AIS, nodes, and unmyelinated axon have elevated sodium conductance
    neuron.mech{t}.iseg.nax = struct(...
        'gbar', 15, ...
        'ena', 50);
    neuron.mech{t}.node.nax = struct(...
        'gbar', 15, ...
        'ena', 50);
    neuron.mech{t}.axon.nax = struct(...
        'gbar', 15, ...
        'ena', 50);
   neuron.mech{t}.all.xtra = struct();
   neuron.mech{t}.all.extracellular = struct();


end

%% Set up cells and run basic simulation with no inputs
cells                    = tree;             % tree morphologies without the source stimulation cells
N                        = 1;                           % Number of Simulations
regions                  = {'basal','proxAp','middleAp','distalAp','tuft'};   % Regions where synapses should be placed
for t                    = 1 : numel (tree)
    % array with as many zeros as nodes in the tree:
    dend{t}              = zeros (size (tree{t}.X), 'double');
    for sim              = 1 : numel (regions)
        % Get ones in the regions you want to activate:
        dend{t}          = dend{t} + double (...
            tree{t}.R (:) == find (strcmp (regions{sim}, tree{t}.rnames)));
    end
    % Get ones in the regions you do not want to activate:
    rest{t}              = double (abs (dend{t}-ones (size (tree{t}.X), 'double')));
    % Calculate the length of each group of regions in order to be able to
    % activate per density:
    dendlength (t)       = sum (len_tree (tree{t}).*dend{t});                    
    restlength (t)       = sum (len_tree (tree{t}).*rest{t});
end

%% Insert netstim for synaptic activity. Set for 3Hz Poisson activity as default
tree{2} = struct('artificial','NetStim','start',10,'interval',333,'number',1000, 'noise', 0.5); %Random synaptic input
tree{3} = struct('artificial','NetStim','start',10,'interval',3,'number',1000, 'noise', 0);%Synchronous synaptic input
t2n_writeTrees (tree,tname, fullfile (treepath, treeFilename));


for t = 1:numel (cells)
   neuronn{t} = neuron;
    iseg_nodes{t} = find(cells{t}.R == find(strcmp(cells{t}.rnames, 'iseg')));
    recnodes = 1;
    neuronn{t}.record{t}.cell = struct('node',recnodes,'record',{'v'});
    nneuron{t}.custom{t}       = [];
end

%% Add simple synapse and connect to netstim

% Find the node corresponding to our chosen distance for the synapse
Pvec = Pvec_tree(tree{t});
dist_vec = Pvec - syn_distance;
for node = 1:length(dist_vec)
    if ~((tree{t}.R(node) == 3 )|| (tree{t}.R(node) == 4) || (tree{t}.R(node) == 5))  
        dist_vec(node) = 99999;
    end
    if dist_vec(node)<=0
        dist_vec(node) = 99998;
    end
end
[~, syn_target] = min(dist_vec);

neuronn{1}.pp{1}.Exp2Syn = struct('node',syn_target,'tau1',0.2,'tau2',2.5,'e',0); %add an excitatory synapse

neuronn{1}.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',syn_target),'delay',0,'threshold',0.5,'weight',0.005);
neuronn{1}.con(2) = struct('source',struct('cell',3),'target',struct('cell',1,'pp','Exp2Syn','node',syn_target),'delay',0,'threshold',0.5,'weight',0.005);

%%

%Files are copied before running t2n to save c and o files of compiled mods
copyfile('./Jarsky_files/lib_mech/', strcat('../Models/', name, '/lib_mech/'), 'f');
copyfile('./Jarsky_files/lib_mech/', './lib_mech/', 'f');

%Make these folders now so that T2N won't get upset
if ~exist('lib_custom', 'dir')
    mkdir('lib_custom');
end
if ~exist('lib_genroutines', 'dir')
    mkdir('lib_genroutines')
end

%This will generate a segmentation fault error; ignore it and save outputs
try
out              = t2n (neuronn,tree, '-d-w-q-m');
catch
end

if ~exist(strcat('../Models/', name, '/Code/sim1/init_cells.hoc'), 'file')
    errordlg(['Model generation for ' name ' failed!']);
    rmpath('Jarsky_files');
    return
end

h = findall(0,'Type','figure','Name','Error in NEURON'); % it returns all the handles for dialog boxes with the title "Error in NEURON"
if isempty(h) % check if such error dialog exists
else
    close(h) % close the error dialog
    disp(['Model generation for ' name ' completed!'])
    disp('--------------------------------------');
end

fclose('all');
%copy lib_mech back to generator to get our c and o files back 
copyfile(strcat('../Models/', name, '/lib_mech/'), './lib_mech/', 'f');



%% Copy across necessary files to model folder
copyfile('./lib_custom/', strcat('../Models/', name, '/lib_custom/'), 'f');
copyfile('./lib_genroutines/', strcat('../Models/', name, '/lib_genroutines/'), 'f');
copyfile('./morphos/', strcat('../Models/', name, '/morphos/'), 'f');

movefile(strcat('../Models/', name, '/Code/sim1/'), strcat('../Models/', name, '/Code/NEURON/'));
delete(strcat('../Models/', name, '/Code/NEURON/neuron_runthis.hoc'));
if exist(strcat('../Models/', name, '/Code/NEURON/tvec.dat'), 'file')
    delete(strcat('../Models/', name, '/Code/NEURON/tvec.dat'));
end

for cell_num = 1:numel(trees)
    if numel(trees) == 1
        copyfile('./TMS_package/', strcat('../Models/', name, '/Code/'), 'f');
    else
        copyfile('./TMS_package/', strcat('../Models/', trees{cell_num}.name, '/Code/'), 'f');
    end
end


%Remove temp folders from the generator folder
delete('./morphos/hocs/*');
rmdir('./morphos/hocs/');
delete('./lib_custom/*');
rmdir('lib_custom');
delete('./lib_genroutines/*');
rmdir('./lib_genroutines/');
delete('./lib_mech/*');
rmdir('./lib_mech/');
rmpath('./Jarsky_files/');
end
