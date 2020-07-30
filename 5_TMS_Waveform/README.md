# TMS Waveform
In this step the waveform of the elcteric field induced by the TMS pulse train is generated. The output waveform file is used later in the NEURON simulations.

## Instructions
Make sure the 'TMS_Waveform' is located within the root NEURON folder, to ensure the output files are located where NEURON code is looking for. Then, Run 'TMS_Waveform.m' in Matlab. A GUI pops up asking for the TMS pulse type (monophasic, biphasic) and then the inter-pulse interval and the number of pulses. After the generation of the pulse train, the resulting waveform is saved in the 'TMS_Waveform_out' directory. That directory includes the following files:
1. 'TMS_E_train.txt' which represent the values of the E-fields waveform (and therefore the quasipotentials due to linear relation) for every sample
2. 'TMS_t_train.txt' which holds the time points of the samples in ms

**Note:** The waveforms are normalized and therefore unitless. The amplitude of the TMS pulse is taken care of during calculation of the spatial distribution of the E-field and quasipotentials.

**For advanced users:** You can create any custom waveform you are interested in by creating the files above with the similar format. Make sure you use the same time steps in generating the waveforms as in the NEURON simulations (default: 0.025 ms).

## Software Requirements
[Matlab](https://www.mathworks.com/) (Tested on Matlab 2019a and 2019b) 
