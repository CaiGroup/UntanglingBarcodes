# -*- coding: utf-8 -*-
"""
Created on Fri Apr  2 09:57:34 2021

@author: jonat
"""

import numpy as np
import pandas as pd
from skimage.morphology import dilation, ball
from skimage.restoration import rolling_ball
from skimage.filters import median
import skimage.io
import tifffile
from scipy.ndimage import shift
import util


hyb = int(snakemake.wildcards['hyb'])

#hyb_im = skimage.io.imread(snakemake.input[0])

offsets = pd.read_csv(snakemake.input[2])
labeled_img = skimage.io.imread(snakemake.input[3])

hyb_stack, md = util.pil_imread(snakemake.input[0], metadata=True)

ind_2_ch_dict = {}
for im_md in md: ind_2_ch_dict[im_md["Channel"]] = im_md["ChannelIndex"] 
#{'640': 0, '561': 1, '488': 2, '405': 3}

ch = int(snakemake.wildcards['ch'])
ch_stack = hyb_stack[ch, :,:,:]
#ch_stack = hyb_stack[int(ind_2_ch_dict[ch]), :,:,:]

#chstack_reshaped = np.zeros([2048,2048,17], np.uint16)
#bgstack_reshaped = np.zeros([2048,2048,17], np.uint16)

bg = util.pil_imread(snakemake.input[1])
#bg = bg[int(ind_2_ch_dict[ch]), :,:,:]
bg = bg[ch, :,:,:]

depth, width, height = np.shape(ch_stack)


#for z in range(depth):
#    chstack_reshaped[:,:,z] = ch_stack[z,:,:]
#    bgstack_reshaped[:,:,z] = bg[z,:,:]

offsets.set_index("hyb", inplace=True)

rbr = snakemake.params['r_ball']
bead_thresh = snakemake.params["bead_thresh"]
bg_sub_multiplier = snakemake.params["bg_sub_multiplier"]

x_t = offsets.loc[hyb, "x"]
y_t = offsets.loc[hyb, "y"]
z_t = offsets.loc[hyb, "y"]

#bgt = shift(bgstack_reshaped, (y_t, x_t, 0), order = 1)
#im_bgt_sub = chstack_reshaped - bg_sub_multiplier*bgt
#im_bgt_sub[chstack_reshaped < bg_sub_multiplier*bgt] = 0

bgt = shift(bg, (0, y_t, x_t), order = 1)
im_bgt_sub = ch_stack - bg_sub_multiplier*bgt
im_bgt_sub[ch_stack < bg_sub_multiplier*bgt] = 0

bead_mask = np.greater(bgt, bead_thresh)
bead_mask = skimage.morphology.binary_dilation(bead_mask, ball(2)) #snakemake.params["r_med_filt"]))

mfiltd = median(im_bgt_sub, ball(snakemake.params["r_med_filt"]))
rb_background = rolling_ball(mfiltd, radius = rbr)#snakemake.params["r_ball"])
im_bgtrb_sub = im_bgt_sub - rb_background
#im_bgtrb_sub[im_bgtrb_sub < 0] = 0
im_bgtrb_sub[im_bgt_sub < rb_background] = 0

im_bgtrb_sub[bead_mask] = 0


mask = labeled_img > 0
print(np.shape(mask))
print(np.shape(im_bgtrb_sub))
for z in range(depth):
    im_bgtrb_sub[z,:,:][mask == 0] = 0

#masked_im_bgtrb_sub = im_bgtrb_sub * mask

#tifffile.imsave(snakemake.output[0], np.int16(im_bgtrb_sub))
skimage.io.imsave(snakemake.output[0], np.int16(im_bgtrb_sub[4,:,:]))

#np.savetxt(snakemake.output[0], im_bgtrb_sub)
