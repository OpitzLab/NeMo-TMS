%This script will run both initialization scripts for T2N and TREES Toolbox
%if they are currently not found to have been installed. Run before
%attempting any model generation in NeMo-TMS.


if ~exist('load_tree', 'file')
    cd ./treestoolbox-master;
    start_trees();
    cd ..
end
if ~exist('t2n_writeTrees', 'file')
    cd ./T2N-master;
    t2n_runthisAfterUnzip();
    cd ..;
end