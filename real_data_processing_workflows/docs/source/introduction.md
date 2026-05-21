# Introduction

This is the documentation for the Snakemake workflows used to process original [SeqFISH+](https://doi.org/10.1038/s41586-019-1049-y) data and a new experiment encoded using Reed-Solomon codes. The workflows use a dots-first approach to decoding with alignment by [fiducial marker matching](https://github.com/CaiGroup/seqfish_fmalign), background subtraction, [ADCG](https://github.com/CaiGroup/SeqFISH_ADCG.jl) dot detection, and [syndrome decoding](https://github.com/CaiGroup/SeqFISHSyndromeDecoding). Both workflows are similar, but with some differences discussed in the following sections and in the manuscript.

Raw image data for the original seqFISH+ experiment is provided on zenodo. Raw image stacks for the Reed-Solomon encoded experiment are too large for a zenodo repository, so compressed preprocessed images are included in the results folder of the Reed-Solomon encoded processing workflow git repository.

For instructions, notebooks, and scripts to produce Reed-Solomon codebooks for seqFISH experiments, go [here](https://github.com/CaiGroup/DisentanglingBarcodesInOpticallyDenseSeqFISHData/tree/main/codebook_generation).

Set parameters for background subtraction and [ADCG](https://github.com/CaiGroup/SeqFISH_ADCG.jl) in the file: [config/config.yaml]

Note: you may get slightly different results when running these workflows when we did. We set the integer linear programming solver to return the best solution it finds after 5 minutes of searching for solutions in large conflicting barcode candidate networks. In these cases, the quality of the approximate solutions found may vary slightly depending on the machine or the version of the optimization solver. We have found variations of around 1/1000.