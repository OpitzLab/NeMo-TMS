% PL_TREE   Topological path length.
% (trees package)
% 
% PL = PL_tree (intree, options)
% ------------------------------
% 
% Returns the topological (!!) path length PL to the root node in the tree.
% For metric path lengths use Pvec_tree.
%
% Input
% -----
% - intree   ::integer:      index of tree in trees or structured tree
% - options  ::string:
%     '-s'   : show
%     {DEFAULT: ''}
%
% Output
% ------
% - PL       ::Nx1 vector:   distances from each node to the root (first
%     node) in the tree 
%
% Example
% -------
% PL_tree      (sample_tree, '-s')
%
% See also   BO_tree LO_tree Pvec_tree
% Uses       ver_tree dA
%
% the TREES toolbox: edit, generate, visualise and analyse neuronal trees
% Copyright (C) 2009 - 2016  Hermann Cuntz

function  PL = PL_tree (intree, options)

% trees : contains the tree structures in the trees package
global       trees

if (nargin < 1) || isempty (intree)
    % {DEFAULT tree: last tree in trees cell array}
    intree   = length (trees);
end

ver_tree     (intree);                 % verify that input is a tree

% use only directed adjacency for this function
if ~isstruct (intree)
    dA       = trees{intree}.dA;
else
    dA       = intree.dA;
end

if (nargin < 2) || isempty (options)
    % {DEFAULT: no option}
    options  = '';
end

% calculating weighted path length:
counter      = 1;
PL           = dA (:, 1);
resPL        = PL;
while sum (resPL == 1) ~= 0
    counter  = counter + 1;
    % use adjacency matrix to walk through tree:
    resPL    = dA * resPL;
    PL       = PL + counter .* resPL;
end
PL           = full (PL);

if strfind       (options, '-s')           % show option
    clf; hold on; 
    HP           = plot_tree (intree, PL, [], [], [], '-b');
    set          (HP, ...
        'edgecolor',           'none');
    colorbar;
    title        ('topological path length');
    xlabel       ('x [\mum]');
    ylabel       ('y [\mum]');
    zlabel       ('z [\mum]');
    view         (2);
    grid         on;
    axis         image;
end

% % shorter but slower (concatenation issue):
% ipar = ipar_tree(index);
% PL = ipar>0;
% PL = sum(PL')-1;




