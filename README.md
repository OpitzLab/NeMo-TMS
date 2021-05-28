# NeMo-TMS Toolbox
Neuron Modeling for TMS (NeMo-TMS) toolbox is for multi-scale modeling of the effects of Transcranial Magnetic Stimulation on single-neuron activity. The instructions and open-source codes provided here enable users with different levels of expertise to investigate the neuronal behavior under TMS.

Refer to the following article for more details:
Shirinpour et al., 2020. Multi-scale Modeling Toolbox for Single Neuron and Subcellular Activity under (repetitive) Transcranial Magnetic Stimulation. bioRxiv, https://doi.org/10.1101/2020.09.23.310219

## Overview
This paradigm incorporates modeling at three scales:

- TMS-induced electric field calculation at the macroscopic scale 
- Simulation of neuronal activity under the external electric field from the previous step
- Simulation of subcellular calcium concentrations based on the membrane voltages calculated in the previous step.

The procedures for running simulations at each scale and the intermediate steps are given in the instructions. This pipeline has been tested on **Windows 10**, and **Linux** (Ubuntu 16.04/18.04). We have tested all the steps except the model generation (step 1) on **macOS Catalina**. However, the codes are developed to run on the native operating system's configurations and possibly not work if there are any modifications (for example, using other forms of shell).

## Instructions
Refer to **Full_Tutorial.pdf** for a comprehensive tutorial on how to use NeMo-TMS.

For a quick guide on how to use NeMo-TMS, refer to **Quick_Guide.pdf**, or open **Quick_Guide.mlapp** in MATLAB (2019b or later).

## Older versions
You can access the older versions of the toolbox from the [Releases](https://github.com/OpitzLab/NeMo-TMS/releases) section.

## Reference
Please cite the following reference:

Shirinpour et al., 2020. Multi-scale Modeling Toolbox for Single Neuron and Subcellular Activity under (repetitive) Transcranial Magnetic Stimulation. bioRxiv, https://doi.org/10.1101/2020.09.23.310219

## Correspondence
aopitz[at]umn.edu (Prof. Alexander Opitz)

