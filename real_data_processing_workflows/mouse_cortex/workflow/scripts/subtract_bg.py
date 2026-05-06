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
from scipy.ndimage import shift

z = int(snakemake.wildcards['z'])
hyb_im = skimage.io.imread(snakemake.input[0])
bg = skimage.io.imread(snakemake.input[1])

offsets = pd.read_csv(snakemake.input[2])
labeled_img = skimage.io.imread(snakemake.input[3])
labeled_img = skimage.transform.resize(labeled_img[:,:], (2048, 2048), order=0)

offsets.set_index("hyb", inplace=True)

hyb = int(snakemake.wildcards['hyb'])
rbr = snakemake.params['r_ball']
bead_thresh = snakemake.params["bead_thresh"]
bead_dilation_radius = snakemake.params["bead_dilation_radius"]


x_t = offsets.loc[hyb, "x"]
y_t = offsets.loc[hyb, "y"]

bgt = shift(bg, (y_t, x_t), order = 1)
im_bgt_sub = hyb_im - bgt
im_bgt_sub[hyb_im < bgt] = 0

bead_mask = np.greater(bgt, bead_thresh)
bead_mask = skimage.morphology.binary_dilation(bead_mask, disk(bead_dilation_radius)) 

mfiltd = median(im_bgt_sub, disk(snakemake.params["r_med_filt"]))
rb_background = rolling_ball(mfiltd, radius = rbr)
im_bgtrb_sub = im_bgt_sub - rb_background
im_bgtrb_sub[im_bgt_sub < rb_background] = 0

im_bgtrb_sub[bead_mask] = 0


mask = labeled_img > 0
masked_im_bgtrb_sub = im_bgtrb_sub * mask

skimage.io.imsave(snakemake.output[0], masked_im_bgtrb_sub)