# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 10:30:17 2021

@author: jonat
"""

import pandas as pd

#sumstats = pd.read_csv(snakemake.input[0])

#for infile in snakemake.input[1:]:
#    sumstats = sumstats.append(pd.read_csv(infile))

#sumstats.to_csv(snakemake.output[0], index=False)

pd.concat([pd.read_csv(infile) for infile in snakemake.input]).to_csv(snakemake.output[0])

