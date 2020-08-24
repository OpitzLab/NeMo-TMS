# Neuron Reconstructions

You can find 10 sample morphology files here. You can use other morphologies as you see fit. However, make sure they have high quality and are suitable for use with the scope of this pipeline.


# Modifications made

The raw ASCII files here have been subjected to some modifications to ensure compatibility with this pipeline.

Thes cells were firstly subjected to a quadratic diameter taper to bring them into concordance with typical diameters for CA1 pyramidal cells. This is due to fluorescence halo effects masking the cells' true diameter. Likewise, the axon diameter is not reliable, and when the artificial myelination algorithm in the model generation step is ran, the axon will be set to a constant diameter.

Additionally, some manual modifications were made to remove any sharp, abrupt kinks or loops in neurite structure, as these are likely reconstruction artefacts.

Finally, a Laplacian smooth (alpha 0.25, 20 iterations) was applied to all regions except for the soma using ProMesh, as this eliminates discontinuities which can cause unrealistic local induced electric field values.
