# -*- coding: utf-8 -*-
"""
Created on Mon May 24 12:45:00 2021

@author: jonat
"""

from scipy.io import loadmat

import numpy as np
import pandas as pd

genes = loadmat(snakemake.input[0])['allNames']
locs1 = loadmat(snakemake.input[1])['tot']
locs2 = loadmat(snakemake.input[2])['tot']
genes = [gene[0][0] for gene in genes]

def reformat(locs):
    gdfs = [] #pd.DataFrame(columns=["x", "y", "fov", "cell", "gene", "gene_num"])
    nfovs, ncells, ngenes = np.shape(locs)
    for fov in range(nfovs):
        for cell in range(ncells):
            for gene_num in range(ngenes):
                if np.size(locs[fov,cell,gene_num]) > 0:
                    fcg_df = pd.DataFrame()
                    fcg_df.loc[:,"x"] = locs[fov,cell,gene_num][:,0]
                    fcg_df.loc[:,"y"] = locs[fov,cell,gene_num][:,1]
                    fcg_df.loc[:,"pos"] = fov
                    fcg_df.loc[:,"cellid"] = cell
                    fcg_df.loc[:,"gene"] = genes[gene_num]
                    fcg_df.loc[:,"gene_num"] = gene_num

                    gdfs.append(fcg_df)

    return pd.concat(gdfs)

cleaned_locs1 = reformat(locs1)
cleaned_locs1.to_csv(snakemake.output[0], index=False)
cleaned_locs2 = reformat(locs2)
cleaned_locs2.to_csv(snakemake.output[1], index=False)
