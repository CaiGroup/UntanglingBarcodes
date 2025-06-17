import skimage.io
import numpy as np
import copy

img = skimage.io.imread(snakemake.input[0])

img_copy = copy.deepcopy(img)
img_copy[1:2046, 1:2046] = 0
edge_regions = np.unique(img_copy)

for reg in edge_regions:
    img[img == reg] = 0

skimage.io.imsave(snakemake.output[0], img)