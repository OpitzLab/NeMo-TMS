# NEURON Simulation
NEURON codes for simulating the neuron behavior under TMS are generated in the previous steps.

## Insturctions
First, generate the desired TMS waveform (refer to TMS_Waveform directory described [here](../5_TMS_Waveform/)). Then ...

For the realistic electric field, place the quasipotentials file calculated [previously](../4_SimNIBS-NEURON-Coupling/) in the root NEURON model folder.

Finally, execute TMS_GUI.hoc, select appropriate parameters, and wait for the simulation to execute.

## Requirements
[NEURON](https://www.neuron.yale.edu/neuron/) (Tested on NEURON 7.5) 
