import skimage
import skimage.io
import tifffile
from scipy.ndimage import shift
import pandas as pd
import numpy as np

labeled_img = skimage.io.imread(snakemake.input[0])

offsets = pd.read_csv(snakemake.input[1])
offsets.set_index("hyb", inplace=True)
hyb = int(snakemake.wildcards["hyb"])
x_t = np.round(offsets.loc[hyb, "x"])
y_t = np.round(offsets.loc[hyb, "y"])

labeled_img = shift(labeled_img, (y_t, x_t), order = 1)


#tifffile.imsave(snakemake.output[0], np.uint8(labeled_img))
skimage.io.imsave(snakemake.output[0], np.uint8(labeled_img))
