# Introduction

Snakemake workflow for processing [SeqFISH+](https://doi.org/10.1038/s41586-019-1049-y) data using alignment by [fiducial marker matching](https://github.com/CaiGroup/seqfish_fmalign), background subtraction, [ADCG] dot detection (as we implement [here](https://github.com/CaiGroup/SeqFISH_ADCG)), and [syndrome decoding] (as we implement [here](https://github.com/CaiGroup/SeqFISHSyndromeDecoding)). 

Set parameters for background subtraction and [ADCG] in the file: [config/config.yaml]