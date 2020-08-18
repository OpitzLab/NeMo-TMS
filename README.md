# Multi-Scale Neuron Modeling
This toolset is for multi-scale modeling of the effects of Transcranial Magnetic Stimulation on single neuron activity. The instructions and open-source codes provided here enable users with a different levels of expertese to investigate the neuronal behavior under TMS.

Refer to the following article for more details:
[JOURNAL INFORMATION]

## Overview
This paradigm incorporates modeling at three scales:

- TMS-induced electric field calculation at the macroscopic scale 
- Simulation of neuronal activity under the external electric field from the previous step
- Simulation of subcellular calcium concentrations based on the membrane voltages calculated in the previous step.

The procedures for running simulations at each scale and the intermediate steps are given in the instructions. This pipeline can run on **Windows**, **Linux**, and **macOS**. However, the codes are developed to run on the native operating system's configurations and possibly not work if there are any modifications (for example, using other forms of shell).

**Note:** throughout this repository, <ins>neuron</ins> refers to the biological neural cells, while <ins>NEURON</ins> refers to the NEURON simulation environment used for the computational modeling of neurons.

## Instructions
To run a full multi-scale model, follow the documentations in each directory in the given order. Make sure to set up the <ins>software requirements</ins> (at the end of the documentation in each folder) before running each step.
1. <code>NEURON-Model-Generation</code> generates neuron models from the morphological reconstructions of neurons (samples in <code>Neuron_Reconstructions</code>)
2. <code>Export_NEURON_Segments</code> exports the coordinates of the NEURON model segments to be used in later steps
3. <code>Electric-Field-Modeling</code> outlines the macroscopic electric field modeling
4. <code>SimNIBS-NEURON-Coupling</code> couples the electric fields in step 3 to the NEURON models generated in step 1
5. <code>TMS_Waveform</code> generates the TMS waveform that is used in the NEURON simulation
6. <code>NEURON-Simulation</code> simulates the activity of the neuron based on the data from previous steps (you can also run this step under the assumption of uniform electric field, in which case, you can skip steps 2 to 4)
7. <code>NEURON_UG4_Interface</code> exports the NEURON simulation data in the compatible format for calcium modeling
8. <code>Calcium-Modeling</code> runs the simulations for calcium concentration based on the data from step 7
9. <code>Visualization</code> outlines the procedure for visualizing the results

## License
[LICENSE STUFF]

Please reference the following article in your publications:
[JOURNAL INFORMATION]

**Correspondence:** aopitz -at- umn.edu (Alex Opitz)
