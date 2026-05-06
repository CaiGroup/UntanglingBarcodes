# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 15:30:08 2021

@author: jonat
"""

import pandas as pd
import numpy as np
import skimage
print(snakemake.input)


offsets = pd.read_csv(snakemake.input[0])

dots = pd.read_csv(snakemake.input[1])

print(dots)

dots['hyb'] = dots['hyb'] + 1
offsets['hyb'] = offsets['hyb'] + 1

dots.set_index("hyb", inplace = True)
offsets.set_index("hyb", inplace = True)


pc_block_table = pd.read_csv(snakemake.input[2])

pc_dict = {row["readout hyb"] : row.pseudocolor for i, row in pc_block_table.iterrows()}
round_dict = {row["readout hyb"] : row["block"] for i, row in pc_block_table.iterrows()}

pcs = []
rounds = []
cellid = []
for i, row in dots.iterrows():
    if i in pc_dict:
        pcs.append(pc_dict[i])
        rounds.append(round_dict[i])
    else:
        pcs.append(None)
        rounds.append(None)

dots["pseudocolor"] = pcs
dots["block"] = rounds

labeled_img = skimage.io.imread(snakemake.input['labeled_image'])
labeled_img = skimage.transform.resize(labeled_img[:,:], (2048, 2048), order=0)

print("np.shape(labeled_img): ", np.shape(labeled_img))
dots = dots.loc[(dots.x < 2048) & (dots.y < 2048) & (dots.x > 0) & (dots.x > 0)]
#dots["cellid"] = [labeled_img[int(row.y)-1, int(row.x)-1] for i, row in dots.iterrows()]
dots["cellid"] = [labeled_img[int(row.y)-1, int(row.x)-1] for i, row in dots.iterrows()]


for hyb in offsets.index:
    #dots.loc[hyb,"x":"y"] -= offsets.loc[hyb,'x_offset':'y_offset']
    if hyb in dots.index and len(dots.loc[hyb]) > 0:
        dots.loc[hyb,"x"] = dots.loc[hyb,"x"] - offsets.loc[hyb,'x']
        dots.loc[hyb,"y"] = dots.loc[hyb,"y"] - offsets.loc[hyb,'y']

dots.reset_index(drop=False, inplace=True)
#dots.sort_values(by=["block", "pseudocolor","x","y"], inplace=True)
dots.to_csv(snakemake.output[0], index=False)
#dots.loc[hyb,"x"] - offsets.loc[hyb,'x_offset']
