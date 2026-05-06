# Simulation

This folder contains the snakemake workflow for our seqFISH+ processing.

The config file specifies:
- The the combinations of number of barcodes and number of random dot for each condition
- the number of replicates
- the width of gaussian dots in simulated images
- the ADCG and decoding parameters

First files containing positions, hybridizations, and identities of dots in each replicate of each condition are generated, then images are generated for these files.
ADCG then attempts to recover dot locations in these files. We use syndrome decoding to decode dots in both the generated dot files and in the ADCG recovered dot locations.
Finally, the decoding results are aggregated.