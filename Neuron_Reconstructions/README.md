# Neuron Reconstruction Files
# Neuron Model Generation
This file will process Neurolucida .ASC files into .SWC files which can be used by the model generation files, and apply a quadratic diameter taper appropriate to a CA1 pyramidal cell.

## Instructions
Before running any scripts utilizing TREES Toolbox, initialize TREES Toolbox by running the included start_trees.m.

First, specify the input .ASC file on line 8 of the taper_unified script.

Next, run the first section of the script only - this is necessary for manual specification of the orientation and the trunk branch point.

Identify what rotation needs to be applied to the tree so that the apical dendrite is pointing downwards in the generated figure. Once this is done, identify the location of the trunk branching point, and insert this point into the c3d variable on line 27 of the script. 

Once both of these manual adjustments is made, the script can be ran in full, and will generate a .SWC file which can be utilized by the model generation phase. The output file will be in accordance with standard SWC practices, with soma as region 1, axon as region 2, basal dendrites as region 3, and apical dendrite as region 4.

## Requirements

MATLAB
TREES Toolbox (http://treestoolbox.org/)
