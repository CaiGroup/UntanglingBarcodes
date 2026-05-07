import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import networkx as nx
import skimage


for infile in snakemake.input:
    results_series_concat_best = pd.read_csv(infile)
    #plt.plot(dotsfirst_stats['est_fp_rate'], dotsfirst_stats['gene_encoding_barcodes'], '.', label="dots first")
    plt.plot(results_series_concat_best.est_fdr, results_series_concat_best['ge'], '.', label = infile)
plt.xlabel('Estimated FDR')
plt.ylabel('Number of decoded gene encoding barcodes')
plt.legend()
plt.tight_layout()
plt.savefig(snakemake.output[0])
plt.close()