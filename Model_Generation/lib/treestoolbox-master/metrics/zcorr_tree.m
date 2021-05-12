% ZCORR_TREE   corrects neurolucida z-artifacts.
% (trees package)
%
% [tree idZ] = zcorr_tree (intree, tZ, options)
% ---------------------------------------------
%
% While reconstructing cells with Neurolucida sudden shifts in the z-axis
% can occur. This function is to correct automatically for those effects.
% Any jump in the z-axis > dz is subtracted from the entire subtree.
%
% Input
% -----
% - intree   ::integer:     index of tree in trees or structured tree
% - tZ       ::value in um: threshold value for counting a wrong dz 
%     {DEFAULT: 5}
% - options  ::string:
%     '-s'   : show
%     '-w'   : waitbar
%     '-m'   : demo movie
%     {DEFAULT: ''}
%
% Output
% -------
% if no output is declared the tree is changed in trees
% - tree     ::structured output tree
% - idZ      ::vector:  output is index of corrected localities.
%
% Example
% -------
% zcorr_tree   (sample_tree, 10, '-s -m') % change nothing since tree good
% % but see:
% zcorr_tree   (sample_tree, 4,  '-s -m')
%
% See also morph_tree flatten_tree
% Uses sub_tree idpar_tree ver_tree
%
% the TREES toolbox: edit, generate, visualise and analyse neuronal trees
% Copyright (C) 2009 - 2016  Hermann Cuntz

function varargout = zcorr_tree (intree, tZ, options)

% trees : contains the tree structures in the trees package
global       trees


if (nargin < 1) || isempty (intree)
    % {DEFAULT tree: last tree in trees cell array}
    intree   = length (trees); 
end;

ver_tree     (intree); % verify that input is a tree structure

% use full tree for this function
if ~isstruct (intree)
    tree     = trees{intree};
else
    tree     = intree;
end

if (nargin < 2) || isempty (tZ)
    tZ       = 5;
end

if (nargin < 3) || isempty (options)
    % {DEFAULT: no option}
    options  = ''; 
end;

ipar             = ipar_tree (tree);
idpar            = ipar (:, 2);
idpar (idpar == 0) = 1;

% compare z to direct parent z:
dZ               = tree.Z (idpar) - tree.Z;
% and index to the nodes at which dZ is bigger than a threshold tZ
idZ              = find (abs (dZ) > tZ);

if strfind       (options, '-m') % show movie option: initialization
    clf;
    shine;
    hold on;
    HP           = plot_tree (tree);
    view         (3);
    grid;
    axis         image;
end

if strfind       (options, '-w') % waitbar option: initialization
    HW           = waitbar (0, 'finding jumps in z node by node...');
    set          (HW, 'Name', '..PLEASE..WAIT..YEAH..');
end

for counter      = 1 : length (idZ)
    if strfind   (options, '-w') % waitbar option: update
        waitbar  (counter / length (idZ), HW);
    end
    isub         = sub_tree (tree, idZ (counter));
    tree.Z (find (isub)) = ...
        tree.Z (find (isub)) + dZ (idZ (counter));
    if strfind   (options, '-m') % show movie option: update
        set      (HP, ...
            'visible',         'off');
        HP       = plot_tree (tree);
        drawnow;
    end
end

if strfind       (options, '-w') % waitbar option: close
    close        (HW);
end

if strfind       (options, '-s') % show option
    clf;
    hold         on;
    HP           = plot_tree (intree);
    set          (HP, ...
        'facealpha',           0.5);
    HP           = plot_tree (tree, [1 0 0]);
    set          (HP, ...
        'facealpha',           0.5);
    HP (1)       = plot (1, 1, 'k-');
    HP (2)       = plot (1, 1, 'r-');
    legend       (HP, ...
        {'before',             'after'});
    set          (HP, ...
        'visible',             'off');
    title        ('z-correction');
    xlabel       ('x [\mum]');
    ylabel       ('y [\mum]'); 
    zlabel       ('z [\mum]');
    view         (3);
    grid         on;
    axis         image;
end

if (nargout > 0) || (isstruct (intree))
    varargout{1}   = tree; % if output is defined then it becomes the tree
else
    trees{intree}  = tree; % otherwise add to end of trees cell array
end

if (nargout > 1)
    varargout{2}   = idZ';
end

