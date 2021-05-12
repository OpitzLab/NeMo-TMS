
[![GitHub tag](https://img.shields.io/github/tag/MarcelBeining/t2n.svg?style=for-the-badge)](https://github.com/MarcelBeining/T2N/releases)
![GitHub top language](https://img.shields.io/github/languages/top/MarcelBeining/t2n.svg?style=for-the-badge)
[![GitHub](https://img.shields.io/github/license/MarcelBeining/t2n.svg?style=for-the-badge)](https://github.com/MarcelBeining/T2N/blob/master/LICENSE)
[![GitHub contributors](https://img.shields.io/github/contributors/MarcelBeining/t2n.svg?style=for-the-badge)](https://github.com/MarcelBeining/T2N/graphs/contributors)
![GitHub repo size in bytes](https://img.shields.io/github/repo-size/MarcelBeining/t2n.svg?style=for-the-badge)
[![GitHub issues](https://img.shields.io/github/issues/MarcelBeining/t2n.svg?style=for-the-badge)](https://github.com/MarcelBeining/T2N/issues)


# Short description
T2N is an extension of the TREES toolbox providing an interface between Matlab and the compartmental modeling environment NEURON.

It was published in (and should be cited this way):  
[Beining M, Mongiat, LA, Schwarzacher SW, Cuntz H, Jedlicka P: T2N as a new tool for robust electrophysiological modeling demonstrated for mature and adult-born dentate granule cells. eLife 2017; 6:e26517](https://elifesciences.org/articles/26517)

- T2N allows an easy generation of real and synthetic morphology single-cell and network models. 
- T2N speeds up simulation time by automatic distributed computing of the simulations (e.g. a series of different current steps, or series of different ion channel blocks) and by using parallel Neuron (e.g. for large-scale networks).
- All mechanisms, point processes (PPs), connections, morphologies and NEURON settings are directly set in a well-defined Matlab structure. 
- For easy location-specific settings and manipulations T2N uses the TREES nodes which are automatically translated into NEURON sections and segments. 
- Matlab and TREES provide convenient ways to analyze the structured simulation output of T2N, thus making T2N a valuable tool for extensive in silico structure-function analyses. 

# More information
For more information on installation and usage, please read the paper and the [Documentation_T2N.docx](https://github.com/MarcelBeining/T2N/blob/master/Documentation%20T2N.docx?raw=true)

# Help
The easiest way is to open an [issue](https://github.com/MarcelBeining/T2N/issues) on github.

# Contributing
I am always open for new ideas and contributions. Please contact me via beining@fias.uni-frankfurt.de

# License
This software is published under the MIT license. For further information read the LICENSE file.
