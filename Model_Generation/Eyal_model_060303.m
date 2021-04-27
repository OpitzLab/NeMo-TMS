%%Eyal and Segev model for TMS - cell0603_03_model_476
% https://elifesciences.org/articles/16553
%  Parameters are specific to cell model!!! 
% Nicholas Hananeia, 2020

function[] = Eyal_model_060303()
% initialize model folder hierarchy in current folder:
addpath('Eyal files');
tstop                    = 1100;%40000;
dt                       = 0.05;

% Define standard parameters
neuron.params            = [];
neuron.params.celsius    = 35;
neuron.params.v_init     = -85;
neuron.params.prerun     = 0;
neuron.params.tstop      = 1500;
neuron.params.dt         = dt;
neuron.params.nseg       = 'd_lambda';
neuron.params.dlambda    = 0.1; neuron.params.freq       = 500;
resampling_interval = 0.5;
syn_distance = 100;


%fileid = strcat('./morphos/', input_cell);
fileid = './morphos/Cell_060303.swc';
[~, name, ~] = fileparts(fileid);
trees = load_tree(fileid); %Specify input morphology here!

axon_type = 5; %menu('Choose desired axon:','Original axon','No axon','Eyal axon','Myelinated stick axon', 'Myelinated original axon');

%%
for cell_num = 1:length(trees)
simName = inputdlg('Enter simulation name; leave blank to use cell filename', 's');
    
    %% Load morphology
treeFilename = './morphos/place_tree.mtr'; %Input file here!
treepath = '';

if length(simName{1})~=0 
    name = simName{1};
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

tree{1} = rot_tree(tree{1}, [0 0 5]);

    
%% Deal with the axon

if axon_type == 2
    tree{1} = strip_axon(tree{1});
elseif axon_type == 3
    tree{1} = addnewaxon(tree{1}, 1.12, 1.12);
elseif axon_type == 4
    tree{1} = strip_axon(tree{1}); %get rid of original axon first
    tree{1} = add_axon(tree{1});    
elseif axon_type == 5
    tree{1} = myelinate_axon(tree{1});
end

if axon_type == 1  
    tree{1}.rnames = {'soma', 'axon', 'basal', 'apical'};
elseif axon_type == 2
    tree{1}.rnames = {'soma', 'axon' 'basal' 'apical'};
elseif axon_type == 3
    tree{1}.rnames = {'soma' 'axon' 'basal' 'apical'};
else   %In this case we need myelin
    tree{1}.rnames = {'soma', 'axon',  'basal', 'apical', 'hill', 'iseg', 'myelin', 'node'};
end

%% For the case of myelinated axon, flatten non-myelin regions to be generic axon
if axon_type == 4 || axon_type == 5
    for i = 1:length(tree{1}.rnames)
        if(strcmp(tree{1}.rnames(i), 'axon'))
            axonR = i;
        elseif strcmp(tree{1}.rnames(i), 'iseg')
            isegR = i;
        elseif strcmp(tree{1}.rnames(i), 'hill')
            hillR = i;
        elseif strcmp(tree{1}.rnames(i), 'node')
            nodeR = i;
        elseif strcmp(tree{1}.rnames(i), 'myelin')
            myelinR = i;
        end
    end
    
    for i  = 1:length(tree{1}.R)
        if tree{1}.R(i) == isegR || tree{1}.R(i) == hillR || tree{1}.R(i) == nodeR
            tree{1}.R(i) = axonR;
        end
        if tree{1}.R(i) == myelinR
            tree{1}.R(i) = 5; %Since we've removed regions, set myelin to be highest
        end
    end
    tree{1}.rnames = {'soma', 'axon', 'basal', 'apical', 'myelin'};
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
            tname = strcat('Eyal_',name);
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
cm                       = 0.488;              % Membrane capacitance (µF/cm²)
Ra                       = 281.78;               % Cytoplasmic resistivity (ohm*cm)
Rm                       = 21406;             % Membrane resistance (ohm/cm²) (uniform)
gpas                     = 1;
e_pas                    = -81.1917;
F                        = 1.9;               %Spine scaling factor
StepDist                 = 60;                %Spine step distance
myelinFactor             = 0.01;               %Myelin capacitance factor
for t                    = 1 : numel (tree)
    neuron.mech{t}.all.pas      = struct ( ...
        'cm',        cm,  ...
        'Ra',        Ra,  ...
        'g',         gpas / Rm,  ...
        'e',         e_pas);
end

% To get the regions that should be ranged:
% taken the regions from tree 1 because all of them have the same region
% definition:
treeregions              = tree{1}.rnames;
noregions                = {'soma', 'axon', 'basal', 'apical'};
x                        = false (size (treeregions));
for r                    = 1 : numel (noregions)
    % <-- Flag the ones that noregions{r} matches:
    x                    = x | ~cellfun (@isempty, strfind (treeregions, noregions{r}));
end
treeregions (x)          = [];    % <-- Delete all the flagged lines at once


%% Scale spines
Pvec = Pvec_tree(tree{1});
Fvec = ones(length(Pvec), 1).*F;
for i = 1:length(Pvec)
    if(Pvec(i) <= StepDist) %Exclude stuff less than 60um from soma
        Fvec(i) = NaN;
        Gvec(i) = NaN;
    end
    if(~(tree{1}.R(i) == 3 || tree{1}.R(i) == 4)) %Only apply to apical/basal dendrites
        Fvec(i) = NaN;
        Gvec(i) = NaN;
    end
end  


%% Add active mechanisms
for t                    = 1 : numel (tree)
    %Spine scaling here
    neuron.mech{t}.range.pas = struct(...
        'g', Fvec./Rm, ...
        'cm', Fvec.*cm);

    neuron.mech{t}.soma.na_ion.ena = 67.5;
    neuron.mech{t}.axon.na_ion.ena = 67.5;
    neuron.mech{t}.soma.k_ion.ek = -102;
    neuron.mech{t}.axon.k_ion.ek = -102;
    

    neuron.mech{t}.soma.NaTg = struct(...
        'gNaTgbar', 0.17392, ...
        'slopem', 13.9283, ...
        'vshiftm', 8.00014, ...
        'vshifth', 11.0455);
    neuron.mech{t}.soma.SK_E2 = struct('gSK_E2bar', 0.0913823);
    neuron.mech{t}.soma.SKv3_1 = struct('gSKv3_1bar', 0.0929495);
    neuron.mech{t}.soma.Ca_LVAst = struct('gCa_LVAstbar', 0.000996713);
    neuron.mech{t}.soma.Ca = struct('gCabar', 0.000702806);
    neuron.mech{t}.soma.CaDynamics_E2 = struct(...
        'gamma', 0.000763324,...
        'decay', 164.168);
    neuron.mech{t}.soma.Nap_Et2 = struct('gNap_Et2bar', 2.13807e-06);
    neuron.mech{t}.soma.K_Pst = struct('gK_Pstbar', 3.66269e-07);
    neuron.mech{t}.soma.K_Tst = struct('gK_Tstbar', 0.0480351);
    neuron.mech{t}.soma.Im = struct('gImbar', 0.000262953);
    
    neuron.mech{t}.axon.NaTg = struct(...
        'gNaTgbar', 4.93565, ...
        'slopem', 14.9916, ...
        'vshiftm', 8.10747, ...
        'vshifth', 0.0102124);
    neuron.mech{t}.axon.SK_E2 = struct('gSK_E2bar', 0.0126023);
    neuron.mech{t}.axon.SKv3_1 = struct('gSKv3_1bar', 1.83291);
    neuron.mech{t}.axon.Ca_LVAst = struct('gCa_LVAstbar', 0.000999427);
    neuron.mech{t}.axon.Ca = struct('gCabar', 0.000985442);
    neuron.mech{t}.axon.CaDynamics_E2 = struct(...
        'gamma', 0.04294,...
        'decay', 466.057);
    neuron.mech{t}.axon.Nap_Et2 = struct('gNap_Et2bar', 6.70607e-05);
    neuron.mech{t}.axon.K_Pst = struct('gK_Pstbar', 0.00196954);
    neuron.mech{t}.axon.K_Tst = struct('gK_Tstbar', 0.000561139);
    neuron.mech{t}.axon.Im = struct('gImbar', 0.000883399);
    
        %Myelin has reduced membrane capacitance, otherwise same properties
        %as dendrites
        if axon_type == 4 || axon_type == 5
         neuron.mech{t}.myelin.pas = struct(...
             'cm', 0.01, ...
             'Ra', 100, ...
             'g', 1/1.125e6);
        end
    

   neuron.mech{t}.all.xtra = struct();
   neuron.mech{t}.all.extracellular = struct();


end

%% Set up cells and run basic simulation with no inputs
cells                    = tree;             % tree morphologies without the source stimulation cells
N                        = 1;                           % Number of Simulations
regions                  = tree{1}.rnames;
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
        dist_vec(node) = 99998;
    end
    if dist_vec(node)<=0
        dist_vec(node) = 99999;
    end
end
[~, syn_target] = min(dist_vec);

neuronn{1}.pp{1}.Exp2Syn = struct('node',syn_target,'tau1',0.2,'tau2',2.5,'e',0); %add an excitatory synapse

neuronn{1}.con(1) = struct('source',struct('cell',2),'target',struct('cell',1,'pp','Exp2Syn','node',syn_target),'delay',0,'threshold',0.5,'weight',0.005);
neuronn{1}.con(2) = struct('source',struct('cell',3),'target',struct('cell',1,'pp','Exp2Syn','node',syn_target),'delay',0,'threshold',0.5,'weight',0.005);

%%

%Files are copied before running t2n to save c and o files of compiled mods
copyfile('./Eyal files/lib_mech/', strcat('../Models/', name, '/lib_mech/'), 'f');
copyfile('./Eyal files/lib_mech/', './lib_mech/', 'f');

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

h = findall(0,'Type','figure','Name','Error in NEURON'); % it returns all the handles for dialog boxes with the title "Error in NEURON"
if isempty(h) % check if such error dialog exists
    errordlg(['Model generation for ' name ' failed!']);
else
    close(h) % close the error dialog
    disp(['Model generation for ' name ' completed!'])
    disp('--------------------------------------');
end

%copy lib_mech back to generator to get our c and o files back 
copyfile(strcat('../Models/', name, '/lib_mech/'), './Eyal files/lib_mech/', 'f');

end

%% Copy across necessary files to model folder
copyfile('./lib_custom/', strcat('../Models/', name, '/lib_custom/'), 'f');
copyfile('./lib_genroutines/', strcat('../Models/', name, '/lib_genroutines/'), 'f');
copyfile('./morphos/', strcat('../Models/', name, '/morphos/'), 'f');

movefile(strcat('../Models/', name, '/Code/sim1/'), strcat('../Models/', name, '/Code/NEURON/'));

for cell_num = 1:numel(trees)
    if numel(trees) == 1
        copyfile('./TMS package/', strcat('../Models/', name, '/Code/'), 'f');
    else
        copyfile('./TMS package/', strcat('../Models/', trees{cell_num}.name, '/Code/'), 'f');
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
rmpath('./Eyal files/');