# Visualization
In this step, we discuss how to visualize the results of the NEURON and calcium concentration simulations. Note that this is just one method for visualizaiton and the users are free to implement their own visualization procdure as they see fit.

## Instruction
run <code>visualize_adjacent()</code> to generate the video of the volatge traces (NEURON simulation results), calcium concentrations, and a side by side representation of both. If you are interested in generating only one of those videos, you can also run <code>visualize_neuron()</code> or <code>visualize_calcium()</code> individually.

## Inputs
<code>visualize_neuron(NEURON_results,output_folder)</code> where:
* **output_folder** points to the folder where the results are generated in
* **NEURON_results** points to the folder including the results of the NEURON simulation. The following files should be present in this folder:
    1. **'voltage_trace.dat'** which contains the voltage data for all segments across time
    2. **'locs\locs_all_seg.txt'** which contains the list of segment coordinates and the connections
    
<code>visualize_calcium(calcium_results,output_folder)</code> where:
* **output_folder** points to the folder where the results are generated in
* **calcium_results** points to the folder including the results of the calcium simulation. The following files should be present in this folder:
    1. **'fullDataOut.txt'** which contains the calcium concentrations for all segments across time
    2. **'outDom.txt'** which contains the list of segment coordinates used in 'fullDataOut.txt'
    3. **neuron_out.swc** is the morphology file used in the calcium simulations.
    
<code>visualize_adjacent(NEURON_results,calcium_results,output_folder)</code> where:
* both the calcium concentration data and voltage traces from NEURON is needed. See above for the requirements.

**Note:** In the above functions, 'output_folder' is created automatically if it doesn't exist.

## Outputs
The following results can be found in the 'output_folder' after running the functions.
* **video_adjacent** contains the video of volatge traces and calcium concentrations together. The top and bottom panels represent the voltage in *mv* and *mol/liter* respectively. 
* **video_neuron** contains the video of volatge traces.
* **video_calcium** contains the video of calcium concentrations.
* **gmsh** folder includes the data used for generating the videos:
    * **png_neuron** contains the snapshots of the voltage distribution over the neuron at each sample.
    * **png_calcium** contains the snapshots of the calcium concentrations over the neuron at each sample.
    * **msh and geo files** are the files that generate the corresponding snapshots in gmsh.

## Software Requirements
[Matlab](https://www.mathworks.com/) (Tested on Matlab 2019b), [Gmsh](https://gmsh.info/) (tested on Gmsh 3.0.6)
