function direction = dir_tree(tree,options)
% -n unit direction vector

if nargin < 2
     options = '';
end

idpar = idpar_tree(tree);

direction =  zeros(numel(tree.X),3);
for n = 1:numel(tree.X)

         direction(n,1) = tree.X(n) - tree.X(idpar(n)); % node to parent node differences
         direction(n,2)  = tree.Y(n) - tree.Y(idpar(n));
         direction(n,3)  = tree.Z(n) - tree.Z(idpar(n));

         if ~isempty(strfind(options,'-n'))
             direction(n,:) = direction(n,:) / norm(direction(n,:));
         end
end