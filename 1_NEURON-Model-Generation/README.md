# Neuron Model Generation
This stage of the pipeline will generate NEURON .hoc model files from .SWC morphology files, or TREES Toolbox generated .MTR files which contain multiple morphologies.

Currently, the pipeline supports the Jarsky model of the CA1 pyramidal cell. 

## Instructions
First, compile the mod files located in Generator/lib_mech/mods. This is a necessary step for each computer that the simulations are executed on. This is done using NEURON's mknrndll. After this is done, move the nrnmech.dll file to Generator/lib_mech.

Next, place the input morphology file (either .SWC or .MTR if multiple morphologies are to be rpocessed) into Generator/morphos alongside the example_tree.mtr file.

Open Jarsky_model.m in MATLAB and change the line labelled "Input file here" to refer to the selected input morphology. Run the script.

If the input morphology lacks an axon, uncomment line 41, and a simple "stick" axon will be added. If the input morphology lacks a soma, uncomment line 40, and a soma will be added at the intersection of the branches.

For each morphology given, a figure will be produced showing the assigned model regions.

Every time the generation is completed, the script will return an error - ignore this: It is expected. After completion, the script will copy all other necessary files to the model folders.

## Requirements

MATLAB
TREES Toolbox (http://treestoolbox.org/)
T2N (https://www.treestoolbox.org/T2N.html)
