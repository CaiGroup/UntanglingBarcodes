# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 11:48:25 2021

@author: jonat
"""

import skimage.io
import numpy as np
import tifffile
import util

stack, md = util.pil_imread(snakemake.input[0], metadata=True)
masks = skimage.io.imread(snakemake.input[1])

ind_2_ch_dict = {}
for im_md in md: ind_2_ch_dict[im_md["Channel"]] = im_md["ChannelIndex"] 
#{'640': 0, '561': 1, '488': 2, '405': 3}
#print(ind_2_ch_dict)

#ch = snakemake.wildcards["ch"]
#ch_stack = stack[int(ind_2_ch_dict[ch]), :,:,:]

ch = int(snakemake.wildcards["ch"])
ch_stack = stack[ch, :,:,:]

depth, width, height = np.shape(ch_stack)


#chstack_reshaped = np.zeros([2048,2048,depth], np.uint16)
#print("mask shape")
#print(np.shape(masks))
#print("ch_stack shape")
#print(np.shape(ch_stack))
#ch_stack[masks > 0] = 0
for z in range(depth):
    ch_stack[z,:,:][masks > 0] = 0

#for z in range(depth):
#    ch_stack[z,:,:][masks[z,:,:] > 0] = 0
#    chstack_reshaped[:,:,z] = ch_stack[z,:,:]

skimage.io.imsave(snakemake.output[0], ch_stack[4,:,:]) 
#tifffile.imwrite(snakemake.output[0], ch_stack) 
#tifffile.imwrite(snakemake.output[0], chstack_reshaped) 
