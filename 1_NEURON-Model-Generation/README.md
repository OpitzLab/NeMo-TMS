# Neuron Model Generation
This stage of the pipeline will generate NEURON HOC model files from .SWC morphology files, or TREES Toolbox generated .MTR files which contain multiple morphologies.

Currently, the pipeline supports the Jarsky model of the CA1 pyramidal cell. (Jarsky T, Roxin A, Kath WL, Spruston N (2005) Conditional dendritic spike propagation following distal synaptic activation of hippocampal CA1 pyramidal neurons. Nat Neurosci 8:1667-76)

## Software Requirements

[MATLAB](https://www.mathworks.com/), [TREES Toolbox](http://treestoolbox.org/), [T2N](https://www.treestoolbox.org/T2N.html), [NEURON](https://www.neuron.yale.edu/neuron/)

Before running any scripts utilizing TREES Toolbox and T2N, initialize them by running <code>start_trees.m</code> and <code>t2n_runthisAfterUnzip.m</code> included in the corresponding folders, otherwise you will receive errors that MATLAB cannot find certain functions. This only needs to be done once per computer.

## Instructions

First, compile the mod files located in 'TMS_Jarksky/Generator/lib_mech/'. This is a necessary step for each computer that the simulations are executed on. Refer to this [link](https://www.neuron.yale.edu/phpBB/viewtopic.php?t=3263) for more information on how to do this.

The morphology files for this script must be in either .SWC or .MTR format. TREES Toolbox contains tools to process Neurolucida .ASC files into these formats.

Given a morphology specified in a file 'cell.asc', the following commands in MATLAB with TREES Toolbox will save it as an '.swc' file:

<code>tree = neurolucida_tree('cell.asc');</code>  
<code>swc_tree(tree, 'cell.swc');</code>

Place the input morphology file (either .SWC or .MTR if multiple morphologies are to be processed) into 'Generator/morphos' alongside the 'place_tree.mtr' file. You can use the sample [morphology files](../Neuron-Reconstructions) provided.

The input morphology must be in standard SWC format, with soma as region 1, axon as region 2, basal dendrites as region 3, and apical dendrite as region 4.

Ensure that MATLAB is open with Generator as the working folder, and run the function <code>Jarsky_model('cell.swc')</code>, where <code>cell.swc</code> can be substituted for your morphology of choice.

A prompt will open asking the user how they wish for the axon to be handled. "Do not alter" will make no changes. "No axon" will strip any existing axon from the morphology. "Stick axon" will add a simple straight axon in the negative Y direction. "Myelinated axon" will apply a basic myelination algorithm with nodes of Ranvier at 100um intervals and at every branch point.

For each morphology given, a figure will be produced showing the assigned model regions.

Click on 'Yes' if you are prompted to initialize model folders.

After completion, the script will generate all necessary files to the 'Model' folder where the 'Generator' directory was located.
