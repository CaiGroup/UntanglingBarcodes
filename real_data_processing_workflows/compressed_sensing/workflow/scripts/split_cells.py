import skimage.io
import pandas as pd
import numpy as np
import copy

img = skimage.io.imread(snakemake.input[0])
bnding_bxs = pd.read_csv(snakemake.input[1])
lbled_im = skimage.io.imread(snakemake.input[2])
offsets = pd.read_csv(snakemake.input[3])
hyb = int(snakemake.wildcards['hyb'])


bnding_bxs.set_index("cell_num", inplace=True)
print(np.shape(img))
print(bnding_bxs)
for outfname in snakemake.output:
    print(outfname)
    cell_dot_diff = outfname.split("_")[-1]
    cell = int(cell_dot_diff.split(".")[0])
    print("cell: ", cell)
    xmin = np.max([0, bnding_bxs.loc[cell, "xmin"] - 1 + int(np.round(offsets.loc[hyb, "y"]))])
    xmax = np.min([bnding_bxs.loc[cell, "xmax"] - 1 + int(np.round(offsets.loc[hyb, "y"])), 2047])
    ymin = np.max([0, bnding_bxs.loc[cell, "ymin"] - 1 + int(np.round(offsets.loc[hyb, "x"]))])
    ymax = np.min([bnding_bxs.loc[cell, "ymax"] - 1 + int(np.round(offsets.loc[hyb, "x"])), 2047])
    img_copy = copy.deepcopy(img)
    img_copy[lbled_im != cell] = 0
    cell_im = img_copy[xmin:xmax, ymin:ymax]
    skimage.io.imsave(outfname, cell_im)



'''
import pandas as pd

fov_dots = pd.read_csv(snakemake.input[0])


for outfname in snakemake.output:
    print(outfname)
    cell_dot_diff = outfname.split("_cell")[1]
    cell = int(cell_dot_diff.split(".")[0])
    print("cell: ", cell)
    cell_dots = fov_dots.loc[fov_dots.cellid == cell]
    cell_dots.to_csv(outfname)
'''