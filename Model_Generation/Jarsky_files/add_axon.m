function tree = add_axon(tree,type)
% Adds an axon to the morphology. Soma region needs to be already defined.

if (nargin < 2) || isempty (type)    % if the number of function input arguments is less than 2, or the input 'type' is empty
    % {DEFAULT: case 'MainenSejnowski'}           
    type  = 'MainenSejnowski';
end

switch type
    case 'MainenSejnowski'
        if ~any(strcmp(tree.rnames,'hill'))
            tree.rnames = [tree.rnames,'hill'];
        end
        if ~any(strcmp(tree.rnames,'iseg'))
            tree.rnames = [tree.rnames,'iseg'];
        end
        if ~any(strcmp(tree.rnames,'myelin'))
            tree.rnames = [tree.rnames,'myelin'];
        end
        if ~any(strcmp(tree.rnames,'node'))
            tree.rnames = [tree.rnames,'node'];
        end
        
        n_axon_seg = 6;     % Change the number of myelin segments 
        su = surf_tree(tree);
        % equiv_diam = sqrt(sum(su(tree.R == find(strcmp(tree.rnames,'soma'))))/(4*pi));
        % equiv_diam = mean(tree.D((ones(numel(tree.X),1) - (tree.R == find(strcmp(tree.rnames,'soma'))))>0));
        equiv_diam = 1;
        dir = dir_tree(tree,'-n');

        dendritic_dir = [0 1 0];%nanmean(dir)/sqrt(sum(nanmean(dir).^2));
        % iseg_diam = equiv_diam/10;
        iseg_diam = equiv_diam/2; % = 0.5 um
        L_hill = 10;
        idpar = 1;
        for l = 1:L_hill
            d = 4*iseg_diam-3*iseg_diam*(l-1)/(L_hill-1);
            [tree,idpar(end+1)] = insert_tree(tree,[0,1,[tree.X(1) tree.Y(1) tree.Z(1)] - l*dendritic_dir,d,idpar(end)],'nix');
        end
        tree.R(idpar(2:end)) = find(strcmp(tree.rnames,'hill')); % give region same as root
        ind = numel(idpar);

        L_iseg = 15;
        for l = 1:L_iseg
            [tree,idpar(end+1)] = insert_tree(tree,[0,1,[tree.X(idpar(ind)) tree.Y(idpar(ind)) tree.Z(idpar(ind))] - l*dendritic_dir,iseg_diam,idpar(end)],'nix');
        end
        tree.R(idpar(ind+1:end)) = find(strcmp(tree.rnames,'iseg')); % give region same as root
        ind = numel(idpar);
        
        L_node = 1;
        L_myelin = 100;
        
        for n = 1:n_axon_seg
            
            for l = 1:L_myelin
                [tree,idpar(end+1)] = insert_tree(tree,[0,1,[tree.X(idpar(ind)) tree.Y(idpar(ind)) tree.Z(idpar(ind))] - l*dendritic_dir,iseg_diam,idpar(end)],'nix');
            end
            tree.R(idpar(ind+1:end)) = find(strcmp(tree.rnames,'myelin')); % give region same as root
            ind = numel(idpar);
            
            if n < n_axon_seg
                for l = 1:L_node
                    [tree,idpar(end+1)] = insert_tree(tree,[0,1,[tree.X(idpar(ind)) tree.Y(idpar(ind)) tree.Z(idpar(ind))] - l*dendritic_dir,iseg_diam*0.75,idpar(end)],'nix');
                end
                tree.R(idpar(ind+1:end)) = find(strcmp(tree.rnames,'node')); % give region same as root
                ind = numel(idpar);
            end
            
        end
end