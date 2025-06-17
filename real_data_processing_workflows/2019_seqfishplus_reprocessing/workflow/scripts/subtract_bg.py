# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 09:57:34 2021

@author: jonat
"""

import numpy as np
import pandas as pd
from skimage.morphology import dilation, disk
from skimage.restoration import rolling_ball
from skimage.filters import median
import skimage.io
import tifffile
from align_subtract_background_img import align_subtract_background_img
from scipy.ndimage import shift


hyb_im = skimage.io.imread(snakemake.input[0])
bg = skimage.io.imread(snakemake.input[1])

offsets = pd.read_csv(snakemake.input[2])
labeled_img = skimage.io.imread(snakemake.input[3])

offsets.set_index("hyb", inplace=True)

hyb = int(snakemake.wildcards['hyb'])
ch = int(snakemake.wildcards['ch'])
rbr = snakemake.params['r_ball']

#im_bgt_sub = align_subtract_background_img(hyb_im, bg, offsets, hyb)
x_t = offsets.loc[hyb, "x"]
y_t = offsets.loc[hyb, "y"]

bgt = shift(bg, (y_t, x_t), order = 1)
im_bgt_sub = hyb_im - bgt
im_bgt_sub[hyb_im < bgt] = 0

mfiltd = median(im_bgt_sub, disk(snakemake.params["r_med_filt"]))
rb_background = rolling_ball(mfiltd, radius = rbr)#snakemake.params["r_ball"])
im_bgtrb_sub = im_bgt_sub - rb_background
#im_bgtrb_sub[im_bgtrb_sub < 0] = 0
im_bgtrb_sub[im_bgt_sub < rb_background] = 0

mask = labeled_img > 0
masked_im_bgtrb_sub = im_bgtrb_sub * mask

#tifffile.imsave(snakemake.output[0], masked_im_bgtrb_sub)
skimage.io.imsave(snakemake.output[0], masked_im_bgtrb_sub)

#np.savetxt(snakemake.output[0], im_bgtrb_sub)
