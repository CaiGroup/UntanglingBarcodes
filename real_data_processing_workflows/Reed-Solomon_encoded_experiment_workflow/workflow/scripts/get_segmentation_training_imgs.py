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

ind_2_ch_dict = {}
for im_md in md: ind_2_ch_dict[im_md["Channel"]] = im_md["ChannelIndex"] 
#{'640': 0, '561': 1, '488': 2, '405': 3}
#print(ind_2_ch_dict)

#ch = snakemake.wildcards["ch"]
seg_stack = stack[[ind_2_ch_dict['488'], ind_2_ch_dict['405']], :,:,:]

#clip bead intensities
seg_stack[seg_stack > 1000] = 1000

#chstack_reshaped = np.zeros([2048,2048,depth], np.uint16)


#for z in range(depth):
#    chstack_reshaped[:,:,z] = ch_stack[z,:,:]

tifffile.imwrite(snakemake.output[0], seg_stack[:,0,:,:]) 
tifffile.imwrite(snakemake.output[1], seg_stack[:,8,:,:]) 
tifffile.imwrite(snakemake.output[2], seg_stack[:,16,:,:]) 
#tifffile.imwrite(snakemake.output[0], chstack_reshaped) 
