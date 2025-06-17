# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 11:48:25 2021

@author: jonat
"""

import os
import skimage.io
import numpy as np
from skimage.morphology import ball
from skimage.restoration import rolling_ball
from skimage.filters import median
import tifffile
import util

#stack = skimage.io.imread(snakemake.input[0])
stack, md = util.pil_imread(snakemake.input[0], metadata=True)
masks = skimage.io.imread(snakemake.input[1])

rb_radius = snakemake.params["rb_radius"]
mfilt_radius = snakemake.params["mfilt_radius"]

ind_2_ch_dict = {}
for im_md in md: ind_2_ch_dict[im_md["Channel"]] = im_md["ChannelIndex"] 
#{0: '643', 1: '594', 2: '561', 3: '488'}
print(ind_2_ch_dict)

ch = snakemake.wildcards["ch"]
#labeled_im = skimage.io.imread(snakemake.input[1])
#mask = labeled_im > 0 

#ch_stack = stack[:,:,:,int(ind_2_ch_dict[ch])]
ch_stack = stack[int(ind_2_ch_dict[ch]), :,:,:]

print("size(ch_stack): ", np.shape(ch_stack))

#print("channel stack shape: ", np.shape(ch_stack))

mfiltd = median(ch_stack, ball(mfilt_radius))
rb_background = rolling_ball(mfiltd, radius = rb_radius)

ch_rbbgs = ch_stack - rb_background
ch_rbbgs[ch_stack<rb_background] = 0
#masked_ch_rbbgs = ch_rbbgs*mask

tifffile.imwrite(snakemake.output[0], ch_rbbgs) #masked_ch_rbbgs)
