# TMS Waveform
In this step the waveform of the elcteric field induced by the TMS pulse train is generated. The output waveform file is used later in the NEURON simulations. The length of the NEURON simulation will be the same as the waveform generated in this step.

## Software Requirements
[Matlab](https://www.mathworks.com/) (Tested on Matlab 2019a and 2019b) 

## Instructions
After running the [first step](../1_NEURON-Model-Generation/), go to the **'TMS_Waveform'** directory inside each NEURON model ('Model\cell_name\sim1\TMS_Waveform'). Run <code>TMS_Waveform.m</code> in that location in Matlab. A GUI pops up asking for the TMS pulse type (monophasic, biphasic) and then the inter-pulse interval and the number of pulses. After the generation of the pulse train, the resulting waveform is saved in the 'TMS_Waveform_out' directory. That directory includes the following files:
1. <code>TMS_E_train.txt</code> which represent the values of the E-fields waveform (and therefore the quasipotentials due to linear relation) for every sample
2. <code>TMS_t_train.txt</code> which holds the time points of the samples in ms

**Note:** The waveforms are normalized and therefore unitless. The amplitude of the TMS pulse is taken care of during calculation of the spatial distribution of the E-field and quasipotentials.

**For advanced users:** You can change the delay before and after the pulse train by adjusting the values of **delay_start** and **delay_end** (units are in *ms*) in line 46 and 47. You can also create any custom waveform you are interested in by creating the files above with the same format. Make sure you use the same time steps in generating the waveforms as in the NEURON simulations (default: 0.025 *ms*).
