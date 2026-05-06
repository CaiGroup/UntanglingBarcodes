# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 10:30:17 2021

@author: jonat
"""

import pandas as pd

#sumstats = pd.read_csv(snakemake.input[0])
sumstats = []
for infile in snakemake.input:
    ss = pd.read_csv(infile)
    ss['z'] = infile.split('z_')[1].split('_')[0]
    sumstats.append(ss)

#sumstats.to_csv(snakemake.output[0], index=False)

pd.concat(sumstats).to_csv(snakemake.output[0])

