# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 10:30:17 2021

@author: jonat
"""

import pandas as pd

replicates = [] #pd.read_csv(snakemake.input[0])

for infile in snakemake.input:
    sim_results = pd.read_csv(infile)
    rep = int(infile.split("rep_")[1].split("_")[0])
    sim_results["rep"] = rep
    replicates.append(sim_results)

replicates_cat = pd.concat(replicates) #[pd.read_csv(infile) for infile in snakemake.input])

replicates_cat.to_csv(snakemake.output[0], index=False)
