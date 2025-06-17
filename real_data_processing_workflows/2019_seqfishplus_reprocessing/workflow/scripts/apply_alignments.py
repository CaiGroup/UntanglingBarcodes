# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 15:30:08 2021

@author: jonat
"""

import pandas as pd
import numpy as np


offsets = pd.read_csv(snakemake.input[0])

dots = pd.read_csv(snakemake.input[1])

dots['hyb'] = dots['hyb'] + 1
offsets['hyb'] = offsets['hyb'] + 1

dots.set_index("hyb", inplace = True)

offsets.set_index("hyb", inplace = True)

#dots.drop(0, inplace=True)
#offsets.drop(0, inplace=True)
dots.drop(81, inplace=True)
offsets.drop(81, inplace=True)


for hyb in offsets.index:
    #dots.loc[hyb,"x":"y"] -= offsets.loc[hyb,'x_offset':'y_offset']
    dots.loc[hyb,"x"] = dots.loc[hyb,"x"] - offsets.loc[hyb,'x']
    dots.loc[hyb,"y"] = dots.loc[hyb,"y"] - offsets.loc[hyb,'y']


dots.to_csv(snakemake.output[0])
#dots.loc[hyb,"x"] - offsets.loc[hyb,'x_offset']
