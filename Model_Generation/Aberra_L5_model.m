%%Aberra model for TMS L5 PC

function[] = Aberra_L5_model()
% install T2N and TREES and add path for necessary files
cd ./lib/;
init_t2n_trees();
cd ..;
addpath('./Aberra_files/');


tstop                    = 1100;%40000;
dt                       = 0.05;

% Define standard parameters
neuron.params            = [];
neuron.params.celsius    = 35;
neuron.params.v_init     = -70;
neuron.params.prerun     = 200;
neuron.params.tstop      = tstop;
neuron.params.dt         = dt;
neuron.params.nseg       = 'dlambda';
neuron.params.dlambda    = 0.025;neuron.params.freq       = 500;


%fileid = strcat('./morphos/', input_cell);
fileid = './morphos/Aberra_human_L5.swc';
axon_type = 4;
syn_distance = 10;
[~, name, ~] = fileparts(fileid);
trees = load_tree(fileid); %Specify input morphology here!

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
    disp('Model by that name already exists; exiting...');
    rmpath('Aberra_files');
    return
end

dist_input = inputdlg('Enter distance of synapse from soma on apical dendrite', 'Enter synapse distance');

if length(dist_input{1})~=0
    syn_distance = str2double(dist_input{1});
end


%to de-group the morphologies (if necessary), and have the different tree
%strucutres in the cell array 'tree':
current_dir = strcat('../Models/', name);
t2n_initModelfolders(current_dir);



%Load morphologies. 
if(length(trees) == 1)
    tree{1,1} = resample_tree(trees, 1);
    neuron.params.exchfolder = strcat('../Models/',name, '/Code/');
else
    tree{1,1} = trees{cell_num};
    name = trees{cell_num}.name;
    neuron.params.exchfolder = strcat('../Models/',name, '/Code/');
end

%% divide morphology here
has_myelin = 1;

tree{1}.rnames = {'soma', 'axon', 'basal', 'apical', 'myelin', 'node' 'unmyelin'};

%% Scale diameters to the human equivalent as in Aberra
% axon_factor = 2.453;
% apical_factor = 1.876;
% basal_factor = 1.946;
% somatic_factor = 2.453;
% 
% for index = 1:length(tree{1}.R)
%     if tree{1}.R(index) == 3
%         tree{1}.D(index) = tree{1}.D(index)*basal_factor;
%     elseif tree{1}.R(index) == 4
%         tree{1}.D(index) = tree{1}.D(index)*apical_factor;
%     elseif (tree{1}.R(index) > 4) && (tree{1}.R(index) ~= 7)
%         tree{1}.D(index) = tree{1}.D(index)*axon_factor;
%     end    
% end



%% Convert the tree into NEURON
for t                    = 1 : numel (tree)
    if ~all (cellfun (@(x) isfield (x, 'NID'), tree)) || ...
            ~all (cellfun (@(x) exist (fullfile ( ...
            current_dir, 'morphos', 'hocs', [x.NID, '.hoc']), 'file'), ...
            tree))
        answer = 'OK';
        if strcmp        (answer, 'OK')
            tree{t}      = sort_tree      (tree{t}, '-LO');
            tname = strcat('Aberra_',name);
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

xplore_tree(tree{1}, '-2');

%% Add passive parameters
cm                       = 1;              % Membrane capacitance (µF/cm²)
Ra                       = 100;               % Cytoplasmic resistivity (ohm*cm)
gpas                     = 3e-5;
e_pas                    = -75;

for t                    = 1 : numel (tree)
    neuron.mech{t}.all.pas      = struct ( ...
        'cm',        cm,  ...
        'Ra',        Ra,  ...
        'g',         gpas,  ...
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


Plen = Pvec_tree(tree{1});
vec = (-0.869600 + 2.087000.*exp((Plen-0.000000).*0.003100)).*0.000080; %For hyperpolarization current
isApical = find(strcmp(tree{1}.rnames, 'apical'));
isApical = tree{1}.R == isApical;
vec(~isApical) = NaN;


%% Add active mechanisms


% ********** Na conductance (gNabar)
nainfo.gnode             = 30.0;                % in S/cm2
nainfo.ena = 67.5;
nainfo.region            = treeregions;
% ********** A-type K+ channel proximal (gAKp) and distal (gAKd)
kainfo.ek                = -102;
kainfo.region            = treeregions;





for t                    = 1 : numel (tree)

    
    neuron.mech{t}.apical.na_ion.ena = 50;
    neuron.mech{t}.apical.pas.cm = 2;
    neuron.mech{t}.apical.k_ion.ek = -85;
    neuron.mech{t}.basal.pas.cm = 2;
    neuron.mech{t}.soma.na_ion.ena = 50;
    neuron.mech{t}.soma.k_ion.ek = -85;
    neuron.mech{t}.axon.na_ion.ena = 50;
    neuron.mech{t}.axon.k_ion.ek = -102;

    neuron.mech{t}.range.Ih = struct('gIhbar', vec);
    neuron.mech{t}.basal.Ih = struct('gIhbar', 0.000080);
    
    
    neuron.mech{t}.apical.NaTs2_t = struct('gNaTs2_tbar', 0.026145);
    neuron.mech{t}.apical.SKv3_1 = struct('gSKv3_1bar', 0.004226);
    neuron.mech{t}.apical.Im = struct('gImbar', 0.000143);
    neuron.mech{t}.apical.Ih = struct('gIhbar', 0);
    
    neuron.mech{t}.axon.NaTa_t = struct('gNaTa_tbar', 3.137968);
    neuron.mech{t}.axon.K_Tst = struct('gK_Tstbar', 0.089259);
    neuron.mech{t}.axon.CaDynamics_E2 = struct(...
        'gamma', 0.002910,...
        'decay', 287.198731);
    neuron.mech{t}.axon.Nap_Et2 = struct('gNap_Et2bar', 0.006827);
    neuron.mech{t}.axon.SK_E2 = struct('gSK_E2bar', 0.007104);
    neuron.mech{t}.axon.Ca_HVA = struct('gCa_HVAbar', 0.000990);
    neuron.mech{t}.axon.K_Pst = struct('gK_Pstbar', 0.973538);
    neuron.mech{t}.axon.SKv3_1 = struct('gSKv3_1bar', 1.021945);
    neuron.mech{t}.axon.Ca_LVAst = struct('gCa_LVAstbar', 0.008752);
    
    neuron.mech{t}.soma.CaDynamics_E2 = struct(...
            'gamma',  0.000609, ...
            'decay', 210.485284);
    neuron.mech{t}.soma.SKv3_1 = struct('gSKv3_1bar', 0.303472);
    neuron.mech{t}.soma.SK_E2 = struct('gSK_E2bar', 0.008407);
    neuron.mech{t}.soma.Ca_HVA = struct('gCa_HVAbar', 0.000994);
    neuron.mech{t}.soma.NaTs2_t = struct('gNaTs2_tbar', 0.983955);
    neuron.mech{t}.soma.Ih = struct('gIhbar', 0.000080);
    neuron.mech{t}.soma.Ca_LVAst = struct('gCa_LVAstbar', 0.000333);
    
    
    if axon_type == 4 
            %Nodes
            neuron.mech{t}.node.NaTa_t = struct('gNaTa_tbar', 2*3.137968);
            neuron.mech{t}.node.K_Tst = struct('gK_Tstbar', 0.089259);
            neuron.mech{t}.node.CaDynamics_E2 = struct(...
                'gamma', 0.002910,...
                'decay', 287.198731);
            neuron.mech{t}.node.Nap_Et2 = struct('gNap_Et2bar', 0.006827);
            neuron.mech{t}.node.SK_E2 = struct('gSK_E2bar', 0.007104);
            neuron.mech{t}.node.Ca_HVA = struct('gCa_HVAbar', 0.000990);
            neuron.mech{t}.node.K_Pst = struct('gK_Pstbar', 0.973538);
            neuron.mech{t}.node.SKv3_1 = struct('gSKv3_1bar', 1.021945);
            neuron.mech{t}.node.Ca_LVAst = struct('gCa_LVAstbar', 0.008752);
			%neuron.mech{t}.node.Im = struct('gImbar'); Only in some morph
			%neuron.mech{t}.node.Ca = struct('gCabar'); Only in some
            
            %Myelin
            neuron.mech{t}.myelin.pas = struct(...
                'cm' , 0.02, ...
                'g_pas' , 1/1.125e6);
    
            neuron.mech{t}.unmyelin.NaTa_t = struct('gNaTa_tbar', 3.137968);
            neuron.mech{t}.unmyelin.K_Tst = struct('gK_Tstbar', 0.089259);
            neuron.mech{t}.unmyelin.CaDynamics_E2 = struct(...
                'gamma', 0.002910,...
                'decay', 287.198731);
            neuron.mech{t}.unmyelin.Nap_Et2 = struct('gNap_Et2bar', 0.006827);
            neuron.mech{t}.unmyelin.SK_E2 = struct('gSK_E2bar', 0.007104);
            neuron.mech{t}.unmyelin.Ca_HVA = struct('gCa_HVAbar', 0.000990);
            neuron.mech{t}.unmyelin.K_Pst = struct('gK_Pstbar', 0.973538);
            neuron.mech{t}.unmyelin.SKv3_1 = struct('gSKv3_1bar', 1.021945);
            neuron.mech{t}.unmyelin.Ca_LVAst = struct('gCa_LVAstbar', 0.008752);
            
    end
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
    %dendlength (t)       = sum (len_tree (tree{t}).*dend{t});                    
    %restlength (t)       = sum (len_tree (tree{t}).*rest{t});
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
    if ~(tree{t}.R(node) == 4)  
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
copyfile('./Aberra_files/lib_mech/', strcat('../Models/', name, '/lib_mech/'), 'f');
copyfile('./Aberra_files/lib_mech/', './lib_mech/', 'f');

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
else
    close(h) % close the error dialog
    disp(['Model generation for ' name ' completed!'])
    disp('--------------------------------------');
end

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
rmpath('./Aberra_files/');

end