# Visualization
In this step, we discuss how to visualize the results of the NEURON and calcium concentration simulations. Note that this is just one method for visualizaiton and the users are free to implement their own visualization procdure as they see fit.

## Software Requirements
[MATLAB](https://www.mathworks.com/) (Tested on MATLAB 2019b), [Gmsh](https://gmsh.info/) (tested on Gmsh 3.0.6)

## Instruction
**Note:** This visualization method requires a lot of time and disk space. Therefore it is only suitable for short simulations (set up in [previous step](../5_TMS_Waveform)).

First, make sure you can call **gmsh** from terminal and through MATLAB (run <code>gmsh</code> in the terminal and run <code>system('gmsh')</code> in MATLAB command window). If you can't, it is usually because you haven't set up the environment variables properly for gmsh. Then, Run <code>visualize_adjacent()</code> to generate the video of the volatge traces (NEURON simulation results), calcium concentrations, and a side by side representation of both. If you are interested in generating only one of those videos, you can also run <code>visualize_neuron()</code> or <code>visualize_calcium()</code> individually. While the code is running, Gmsh is called to generate necessary images. Make sure not to click on anything or minimize Gmsh.

If you get an error similar to 'The specified profile is not valid.', you don't have the MPEG-4 codec. To resolve this, you can either install the proper codec or change the format of the video to one of the existing profiles (edit <code>VideoWriter(...)</code> accordingly in the line the error occurs).

## Inputs
<code>visualize_neuron(NEURON_results,output_folder)</code> where:
* **output_folder** points to the folder where the results are generated in
* **NEURON_results** points to the folder including the results of the NEURON simulation. The following files should be present in this folder:
    - **'voltage_trace.dat'** which contains the voltage data for all segments across time
    - **'locs\locs_all_seg.txt'** which contains the list of segment coordinates and their connections
    
<code>visualize_calcium(calcium_results,output_folder)</code> where:
* **output_folder** points to the folder where the results are generated in
* **calcium_results** points to the folder including the results of the calcium simulation. The following files should be present in this folder:
    - **'fullDataOut.txt'** which contains the calcium concentrations for all segments across time
    - **'outDom.txt'** which contains the list of segment coordinates used in 'fullDataOut.txt'
    - **'neuron_out.swc'** is the morphology file used in the calcium simulations.
    
<code>visualize_adjacent(NEURON_results,calcium_results,output_folder)</code> where:
* Both the calcium concentration data and voltage traces from NEURON is needed. See above for the requirements.
* Make sure you run the calcium simulations with the same timesteps and duration as the voltage trace, otherwise you can't visualize them side-by-side using this function.

**Note:** In the above functions, 'output_folder' is created automatically if it doesn't exist.

## Outputs
The following results can be found in the 'output_folder' after running the functions.
* **video_adjacent** contains the video of volatge traces and calcium concentrations together. The top and bottom panels represent the local membrane voltage in *mv* and calcium concentration in *mol/liter* respectively. 
* **video_neuron** contains the video of volatge traces.
* **video_calcium** contains the video of calcium concentrations.
* **gmsh** folder includes the data used for generating the videos (You can delete this folder after the the script is finished to save space):
    * **png_neuron** contains the snapshots of the voltage distribution over the neuron at each sample.
    * **png_calcium** contains the snapshots of the calcium concentrations over the neuron at each sample.
    * **msh and geo files** are the files that generate the corresponding snapshots in gmsh.
