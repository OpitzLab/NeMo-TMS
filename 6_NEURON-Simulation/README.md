# NEURON Simulation
NEURON codes for simulating the neuron behavior under TMS are generated in the previous steps.

## Insturctions

First, generate the desired TMS waveform (refer to TMS_Waveform directory described [here](../5_TMS_Waveform/)). Then ...

For the realistic electric field, place the quasipotentials file calculated [previously](../4_SimNIBS-NEURON-Coupling/) in the root NEURON model folder.

Then, execute GUI_params.hoc and select the desired parameters. If the realistic electric field condition is chosen, be sure to select it here. Otherwise, choose field amplitude and direction using either polar or Cartesian co-ordinates. When this is completed, it will generate a parameters file.

Finally, run TMS_script.hoc and wait for the simulation to complete. It will generate a voltage trace for every segment in results/voltage_trace.dat, and will save index and location files to results/locs. The simulation will exit when it has finished running.

## Software Requirements
[NEURON](https://www.neuron.yale.edu/neuron/) (Tested on NEURON 7.5) 
