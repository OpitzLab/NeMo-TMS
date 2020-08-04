% QUADFUNCDIAM_TREE   Map quadratic and scale proximal diameter tapering to tree.
% (trees package)
% 
% tree = quadfuncdiam_tree (intree,  parameters, fhandle, options, P, ldend)
% -------------------------------------------------------------------
%
% QUADFUNCDIAM_TREE adds to quaddiameter_tree an overlaying function,
% which scales the diameter starting at the distance parameter(4). The function can
% be any function defined by fhandle.
%
% Input
% -----
% - intree::integer:index of tree in trees or structured tree
% - parameters:: vector 1*4:
%                  parameters(1)= slope of quadratic fit {DEFAULT: 0.5}
%                  parameters(2)= slope of scaling function {DEFAULT: 0.5}
%                  parameters(3)= offset {DEFAULT: 0.5}
%                  parameters(4)= distance from the root {DEFAULT: 0}
% - fhandle::scaling function, which scales the diameter starting at 
%    parameter(4)
% - options::string: {DEFAULT ''}
%    '-s' : show
%    '-w' : waitbar
% - P::matrix of three columns: parameters for the quadratic equation in
%    dependence of the root to tip length given in:
% - ldend::vertical vector, same length as P: typical lengths at which P are 
%    given
%
% Output
% ------
% if no output is declared the tree is changed in trees
% - tree:: structured output tree
%
% Example
% -------
% quadfuncdiam_tree (sample_tree, [], [], '-s')
%
% See also
% quaddiameter_tree
% Uses Pvec_tree ipar_tree T_tree ver_tree dA D

function  varargout = quadfuncdiam_tree2 (intree,  parameters, fhandle, options, P, ldend)
% trees : contains the tree structures in the trees package
global trees

if (nargin < 1)||isempty(intree),
    intree = length (trees); % {DEFAULT tree: last tree in trees cell array}
end;

ver_tree (intree); % verify that input is a tree structure

% use full tree for this function
if ~isstruct (intree),
    tree = trees {intree};
else
    tree = intree;
end

if (nargin < 2)||isempty(parameters),
    parameters(1) = 0.5; % slope of quadratic fit
    parameters(2) = 0.5; % slope of scaling function
    parameters(3) = 0.5; % offset
    parameters(4) = 0; %   distance from the root
end

if (nargin < 3)||isempty(fhandle),
    fhandle = @(x) 0;    % if no function is defined, there is no scaling
end

if (nargin <4)||isempty(options),
    options = ''; % {DEFAULT: no option}
end

if (nargin <5)||isempty(P),
    load P % {DEFAULT: parameters calculated for optimal current transfer for
           % branches on their own}
end

if (nargin <6)||isempty(ldend),
    load ldend % {DEFAULT: length values of branches for which P is given
               % quaddiameter_tree uses the P whos ldend is closest to the
               % path length for each path to termination point}
end

N = size (tree.dA, 1); % number of nodes in tree
Plen = Pvec_tree (tree)'; % path length from the root [um]
tree.D = ones (N, 1) .* 0.5; % first set diameter to 0.5 um
% NOTE! I'm not sure about the following line:
ipari  = [(1 : N)' ipar_tree(tree)]; % parent index structure incl. node itself twice
ipariT = ipari (T_tree (tree), :);   % parent index paths but only for termination nodes

if strfind (options, '-w'), % waitbar option: initialization
    HW = waitbar (0, 'calculating quad diameter...');
    set (HW, 'Name', '..PLEASE..WAIT..YEAH..');
end

Ds = zeros (size (ipariT));
for ward = 1 : size (ipariT, 1);
    if strfind (options, '-w'), % waitbar option: update
        if mod (ward, 500) == 0,
            waitbar (ward / size (ipariT, 1), HW);
        end
    end
    iipariT   = ipariT (ward, ipariT (ward, :) ~= 0);
    iipariT   = fliplr (iipariT);
    pathh     = Plen (iipariT);
    [~, i2]   = min ((pathh (end) - ldend).^2); % find which ldend is closest to path length
    quadpathh = polyval (P (i2, :), pathh) .* parameters(1);
    quadpathh = quadpathh.*(pathh<parameters(4)).*(fhandle(-parameters(2) ...
        *(pathh-parameters(4)))+1)+quadpathh.*(pathh>=parameters(4));
    Ds (ward, 1 : length (quadpathh)) = fliplr (quadpathh); % apply the diameters
end

if strfind (options, '-w'), % waitbar option: close
    close (HW);
end

% average the diameters for overloaded nodes (there might be a better way
% to do this than averaging):
for ward = 1 : N,
    iR = ipariT == ward;
%     tree.D (ward) = mean (Ds (iR));
    tree.D (ward) = max (Ds (iR));
end

tree.D = tree.D + parameters(3); % add offset diameter

if strfind (options, '-s'), % show option
    clf; hold on; plot_tree (tree, [0 0 0]);
    title  ('quadratic diameter tapering');
    xlabel ('x [\mum]'); ylabel ('y [\mum]'); zlabel ('z [\mum]');
    view (3); grid on; axis equal;
end

if (nargout == 1)||(isstruct(intree)),
    varargout {1}  = tree; % if output is defined then it becomes the tree
else
    trees {intree} = tree; % otherwise the orginal tree in trees is replaced
end