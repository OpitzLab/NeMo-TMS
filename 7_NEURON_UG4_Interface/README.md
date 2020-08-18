# NEURON-UG4 Interface
In this step we describe how the simulation results from NEURON can be converted into a format compatible with UG4 for calcium simulations. UG4 requires the following data to work:
* Voltage traces across the neuron segments at each time point
* The morphology file (SWC format) corresponding to the neuron segments (Note that this morphology file may be different from the original morphology file you used to generate NEURON models).

## Software Requirements
[MATLAB](https://www.mathworks.com/) (Tested on Matlab 2019b)

## Instructions
The <code>NEURON_UG4_interface</code> directory and its contents are automaticlly placed inside each NEURON model folder ('Model\cell_name\sim1\NEURON_UG4_interface') when the NEURON model is generated in the [first step](../1_NEURON-Model-Generation/). This ensures the files are in the correct location to access the input data. Run <code>export_data()</code> in this directory in Matlab. The following files will be exported in this foler:
1. <code>neuron_out.swc</code> represents the morphology file.
2. Inside the 'voltage_data_calcium' folder, (<code>vm_#####.dat</code>) files represent the voltage data from the NEURON simulation. Each file includes the voltage data of all segments (and their coordinates) at one simulation time point.
