%% CA1 tapering algorithm: requires specification of trunk branch point.
% Cell must be oriented with apical dendrite facing downwards.
% Steffen Platschek original, modified by Nicholas Hananeia, 2020.
% NOTE: RUN FIRST SECTION SEPARATELY FIRST.

clear all;
close all;
tree = neurolucida_tree('cell.ASC'); %Specify input file here!
swc_tree(tree, 'cell.swc');


cell_name_import = [{'cell'}];
cell_name = cell_name_import{1};
extension = '.swc';

ntree = length(cell_name_import);
wtrees = cell(ntree,1);

tree = load_tree([cell_name_import{1} extension]);
tree = rot_tree(tree, [0 0 0]); %be sure to orient tree so that it's facing down
plot_tree(tree); %Identify the trunk location using this plot before continuing!
wtrees{1} = tree;

%% >>>>MANUALLY SPECIFY APICAL TRUNK ENDING LOCATION HERE!!!!!!!!!!!!!<<<<


c3d = [0 0 0];

%% SECTION 1: REMOVE Z AXIS JUMPS
for itree = 1:ntree
    wtree = tree;
 
    figure(5);
    plot_tree(tree);
    
    threshold_slope = 2;
    threshold_length = 4;

    T = T_tree(wtree);
    T_=find(T);
    max = length(T_);
    ward = 1;
    
    
    while ward <= max
        T = T_tree(wtree);
        T_=find(T);
        ipar = ipar_tree(wtree);
        idpar = idpar_tree(wtree);
        Pvec = Pvec_tree(wtree);
        slope = (wtree.Z(idpar)-wtree.Z)./sqrt((wtree.X(idpar)-wtree.X).^2+(wtree.Y(idpar)-wtree.Y).^2);

        pathT = ipar(T_(ward),ipar(T_(ward),:)~=0);
        pathT = flip(pathT);

        path = Pvec(pathT);
        logi_slope_ = find(abs(slope(pathT))>threshold_slope);
        if(~isempty(logi_slope_))
            i1=0;
            while i1<length(logi_slope_)
                i0 = i1+1;
                i1 = find(logi_slope_(i0+1:end)-logi_slope_(i0:end-1)>1,1,'first')+i1;
                if isempty(i1) 
                    i1=length(logi_slope_);
                end
    %             display([i0 i1])
    %             pause(1);
                if(path(i1)-path(i0)>threshold_length);
                    sub = logical(sub_tree(wtree,pathT(logi_slope_(i1))));
                    delete = pathT(logi_slope_(i0:i1));
                    ind = pathT(logi_slope_(i0)-1);
                    Z_shift = wtree.Z(pathT(logi_slope_(i1))) - wtree.Z(ind);
                    wtree.Z(sub) = wtree.Z(sub)-Z_shift;

                    plot3(wtree.X(ind),wtree.Y(ind),wtree.Z(ind),'bo');

                    plot3(wtree.X(delete),wtree.Y(delete),wtree.Z(delete),'ro');

                    wtree = delete_tree(wtree, delete);
                end
            end
        end
        ward = ward + 1;
        max = length(T_);
    end
    figure(1);
    plot_tree(wtree);
    wtree = repair_tree(wtree);
    title('fixed jumps')

    view(51,22);
    % print_export = './';
    % set(gcf,'PaperPositionMode','auto');
    % tprint       ([print_export cell_name{itree} 'fj'], '-tif', [20 16]);


    Y = wtree.Y;
    X = wtree.X;
    smtree = smooth_tree(wtree,[], .9,20);
    smtree.Y = Y;
    smtree.X = X;


     plot_tree(tree);
     plot_tree(smtree,[1 0 0], [-4 0 10]);
     title('before(black) after(red)')
%     view(51,22);
    % print_export = './';
    % set(gcf,'PaperPositionMode','auto');
    % tprint       ([print_export cell_name{itree} 'before_after'], '-tif', [20 16]);
    wtrees{itree} = smtree;
end



%% SECTION 2: QUADRATIC DIAMETER TAPERING
color = [0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1];
cells = 1;

wtrees{1}.rnames{2} = 'axon';
wtrees{1}.rnames{1} = 'apical_dendrite';
wtrees{1}.rnames{3} = 'soma';
wtrees{1}.rnames{4} = 'basal_dendrite';

    figure(2);
    hold on;
    tmp_tree = wtrees{1};
    tmp_tree.Z = tmp_tree.Z.*0;
    [~, v(1)] = min(eucl_tree(tmp_tree,[c3d(1,1:2) 0]));
    plot_tree(wtrees{1},[],[1*200 0 0]);
    plot(wtrees{1}.X(v(1))+1*200,wtrees{1}.Y(v(1)),'ro','MarkerSize',10);
    hold off;

% vlachos fit for trunk!
% Ptr = [0.16 0 0.6 0];
Ptr = [0.127 0 0.4 0];

% vlachos fit for oblique!
% Pobl = [0.085 0 0.35 0];
Pobl = [0.071 0 0.35 0];

% vlachos basal dendrites
Pbas = [0.07 0.05 0.2662 30];

fhandle = @(x) exp(x) - 1;

trunk_trees = cell(cells,1);
oblique_trees = cell(cells,1);
basal_trees = cell(cells,1);
axon_trees = cell(cells,1);

for ward = 1:cells
    tree = wtrees{ward};
    trunk_trees{ward} = quadfuncdiam_tree2(tree,Ptr,fhandle);
    oblique_trees{ward} = quadfuncdiam_tree2(tree,Pobl,fhandle);
    basal_trees{ward} = quadfuncdiam_tree2(tree,Pbas,fhandle);
    axon_trees{ward} = wtrees{ward}; %%we do not want to taper the axon at all
end

node = v;
for ward = 1:cells;
    
    
    tree2 = wtrees{ward};
    
    
    % trunk!
    tree = trunk_trees{ward};
    ipar = ipar_tree(tree);
    logi = ipar(node(ward),ipar (node(ward), :) ~= 0);
    tmp = zeros(length(tree.dA),1);
    tmp(logi) = 1;
    logi_trunk = logical(tmp);
    
    x_vec = tree.X(logi_trunk);
    y_vec = tree.Y(logi_trunk);
    z_vec = tree.Z(logi_trunk);
    d_vec = tree.D(logi_trunk);
    
    for ward2 = 1:sum(logi_trunk)
        [~, i1 ] = min((tree2.X-x_vec(ward2)).^2+(tree2.Y-y_vec(ward2)).^2 + ...
            (tree2.Z-z_vec(ward2)).^2);
        tree2.D(i1) = d_vec(ward2);
%         plot(tree2.X(i1), tree2.Y(i1),'bo');
%         plot(x_vec(ward2)+20, y_vec(ward2),'ro');
    end
    
    % oblique!
    tree = oblique_trees{ward};
    apicalR = find(strcmp(tree.rnames, 'apical_dendrite'));
    logi = ipar(node(ward),ipar (node(ward), :) ~= 0);

    tmp = zeros(length(tree.dA),1);
    tmp(logi) = 1;
    logi_oblique = ~tmp & tree.R == apicalR;
    
    x_vec = tree.X(logi_oblique);
    y_vec = tree.Y(logi_oblique);
    z_vec = tree.Z(logi_oblique);
    d_vec = tree.D(logi_oblique);
    
    for ward2 = 1:sum(logi_oblique)
        [~, i1 ] = min((tree2.X-x_vec(ward2)).^2+(tree2.Y-y_vec(ward2)).^2 + ...
            (tree2.Z-z_vec(ward2)).^2);
        tree2.D(i1) = d_vec(ward2);
%         plot(tree2.X(i1), tree2.Y(i1),'bo');
%         plot(x_vec(ward2)+20, y_vec(ward2),'ro');
    end
    
    % basal!
    tree = basal_trees{ward};
    basalR = find(strcmp(tree.rnames, 'basal_dendrite'));
    logi_basal = tree.R == basalR;
    
    x_vec = tree.X(logi_basal);
    y_vec = tree.Y(logi_basal);
    z_vec = tree.Z(logi_basal);
    d_vec = tree.D(logi_basal);
    
    for ward2 = 1:sum(logi_basal)
        [~, i1 ] = min((tree2.X-x_vec(ward2)).^2+(tree2.Y-y_vec(ward2)).^2 + ...
            (tree2.Z-z_vec(ward2)).^2);
        tree2.D(i1) = d_vec(ward2);
%         plot(tree2.X(i1), tree2.Y(i1),'bo');
%         plot(x_vec(ward2)+20, y_vec(ward2),'ro');
    end
    
    
    %axon + soma!
    tree = axon_trees{ward};
    axonR = find(strcmp(tree.rnames, 'axon'));
    somaR = find(strcmp(tree.rnames, 'soma'));
    logi_axon = tree.R == axonR;
    logi_soma = tree.R == somaR;
    logi_axon = logi_axon | logi_soma; %Since I don't want either tapered
    wtrees{ward} = tree2;
    x_vec = tree.X(logi_axon);
    y_vec = tree.Y(logi_axon);
    z_vec = tree.Z(logi_axon);
    d_vec = tree.D(logi_axon);
    
    for ward2 = 1:sum(logi_axon)
        [~, i1 ] = min((tree2.X-x_vec(ward2)).^2+(tree2.Y-y_vec(ward2)).^2 + ...
            (tree2.Z-z_vec(ward2)).^2);
        tree2.D(i1) = d_vec(ward2);
%         plot(tree2.X(i1), tree2.Y(i1),'bo');
%         plot(x_vec(ward2)+20, y_vec(ward2),'ro');
    end
    
    wtrees{ward} = tree2;
    
end

figure(3);
hold on;
plot_tree(wtrees{ward});
plot_tree(trunk_trees{ward},[1 0 0], [10 0 0]);
hold off;

%% SECTION 3: EXPORT TO STANDARDIZED SWC
wtrees{1} = rot_tree(wtrees{1}, [0 0 180]); %orient with conventional apical-up
figure(4);
xplore_tree(wtrees{1}, '-2');


trees = wtrees{1};
cells = length(trees);


%save tree as standard SWC
for i = 1:length(trees.D)
   %apical dendrite regions
   if (strcmp(trees.rnames(trees.R(i)), 'apical_dendrite'))
       trees.R(i) = 4;
   %basal dendrite regions
   elseif (strcmp(trees.rnames(trees.R(i)), 'basal_dendrite'))
       trees.R(i) = 3;
   %soma
   elseif (strcmp(trees.rnames(trees.R(i)), 'soma'))
       trees.R(i) = 1;
   elseif (strcmp(trees.rnames(trees.R(i)), 'axon'))
       trees.R(i) = 2;
   end
end

trees = tran_tree(trees, [trees.X(1), trees.Y(1), trees.Z(1)].*(-1));

swc_tree(trees, strcat(cell_name, '_standard.swc'));









