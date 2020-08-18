# NEURON Simulation
NEURON codes for simulating the neuron behavior under TMS are generated in the previous steps.

## Instructions

First, generate the desired TMS waveform (refer to TMS_Waveform directory described [here](../5_TMS_Waveform/)). Then, execute **GUI_params.hoc** (<code>nrngui GUI_params.hoc</code> in terminal, double-click on the file on Windows) to select the desired parameters using the GUI. The realistic electric field condition allows you to simulate neuronal activity under the electric fields calculated in the FEM model. In this case, be sure to select the corresponding quasipotentials file from the [previous step](../4_SimNIBS-NEURON-Coupling/). Otherwise, you can skip steps 2-4 and choose to simulate the neuron under a uniform electric field. In this case, select the field amplitude (in _V/m_) and direction using either polar or Cartesian co-ordinates. The Cartesian co-ordinates will be normalized to unit vector, therefore the magnitude of the vector is ignored. After the user has selected the desired parameters, it will generate a parameters file. 

Finally, run **TMS_script.hoc** (<code>nrniv TMS_script.hoc</code> in terminal, double-click on the file on Windows) to run the NEURON simulation based on the parameters you selected.

**For advanced users:** You can make new parameters files or edit existing ones yourself, as long as it has the same name (params.txt) and you follow the same format that **GUI_params.hoc** generates. Note that the order of the parameters defined in the file matters.

## Outputs
The following output files can be found in the 'results' folder after running the simulation:

<code>tvec.dat</code> saves the time value at every time step in the simulation.

<code>voltage_trace.dat</code> saves the membrane voltage at every segment in the morphology. It is represented as a matrix where every row represents a simulation time step, and every column represents a segment, with the column index being the segment number.

Within the subfolder 'locs':

<code>index.dat</code> saves an index where each segment is associated with its segment and its morphological region. The first column states the region and section identifier, the second column states the diameter, the third states at what point along the section the segment is located, and the final column is the segment number.

<code>locs_all_seg.txt</code> saves the location of each segment as well as the location of its parent segment. Each row refers to a segment, the first 3 columns are the (X,Y,Z) co-ordinates of the segment itself, and columns 4,5, and 6 refer to the (X,Y,Z) co-ordinates of the parent.

<code>diam_all.txt</code> saves the diameter of every section.

## Software Requirements
[NEURON](https://www.neuron.yale.edu/neuron/) (Tested on NEURON 7.5) 
