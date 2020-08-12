# NEURON Simulation
NEURON codes for simulating the neuron behavior under TMS are generated in the previous steps.

## Instructions

First, generate the desired TMS waveform (refer to TMS_Waveform directory described [here](../5_TMS_Waveform/)). Then ...


Then, execute GUI_params.hoc (<code>nrngui GUI_params.hoc</code> in terminal) and select the desired parameters. If the realistic electric field condition is chosen, be sure to select the corresponding quasipotentials file calculated [previously](../4_SimNIBS-NEURON-Coupling/). Otherwise, choose field amplitude and direction using either polar or Cartesian co-ordinates. When this is completed, it will generate a parameters file.


Finally, run TMS_script.hoc (<code>nrniv TMS_script.hoc</code> in terminal) and wait for the simulation to complete.

## Inputs

<code>TMS_sim()</code> is the primary simulation function. It will execute the simulation and return 1 if the cell generated at least one spike, and 0 otherwise.

<code>GUI()</code> will launch a graphical user interface for choosing simulation parameters. It will generate a file <code>params.txt</code> which is necessary for the simulation to run.

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
