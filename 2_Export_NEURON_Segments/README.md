# Export NEURON Segments
This script is automatically placed in the NEURON model directory after running the [first step](../1_NEURON-Model-Generation). The coordinates and other information about the segments implemented in the HOC files are exported from the generated NEURON model. These data are used in multiple subsequent steps.

## Instructions
Run <code>save_locations.hoc</code> (<code>nrniv save_locations.hoc</code> in terminal) in the NEURON model directory. The output data will be exported to the **'results\locs\'** folder.

## Software Requirements
[NEURON](https://www.neuron.yale.edu/neuron/) (Tested on NEURON 7.5) 
