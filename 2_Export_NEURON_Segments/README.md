# Export NEURON Segments
In this step, the coordinates and other information about the segments implemented in the HOC files are exported from the generated NEURON model. These data are used in multiple subsequent steps.

## Software Requirements
[NEURON](https://www.neuron.yale.edu/neuron/) (Tested on NEURON 7.5) 

## Instructions
After running the [first step](../1_NEURON-Model-Generation), the NEURON codes corresponding to each cell morphology are generated inside **'1_NEURON-Model-Generation\TMS_Jarsky\Model\cell_name\sim1'**. In that directory, run <code>save_locations.hoc</code> (double-click on the file on windows, <code>nrniv save_locations.hoc</code> in terminal in Linux and macOS) in the NEURON model directory. The output data will be exported to the **'results\locs\'** folder.
