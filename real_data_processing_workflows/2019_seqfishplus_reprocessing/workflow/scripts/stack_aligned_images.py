# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 10:11:46 2021

@author: jonat
"""

import skimage.io
import numpy as np
import pandas as pd
import tifffile

def translate_img(img, offsets, hyb):
    x_t = offsets.loc[hyb, "x"]
    y_t = offsets.loc[hyb, "y"]

    x_p = x_t % 1
    y_p = y_t % 1

    x_tf = int(np.floor(x_t))
    x_tc = int(np.ceil(x_t))
    y_tf = int(np.floor(y_t))
    y_tc = int(np.ceil(y_t))

    xsize, ysize = np.shape(img)
    translated = np.zeros((xsize, ysize))

    for i in range(xsize):
        for j in range(ysize):
            if j - x_tc < xsize and i - y_tc < ysize and j - x_tc > 0 and i - y_tc > 0 and j - x_tf < xsize and i - y_tf < ysize and j - x_tf > 0 and i - y_tf > 0:
                translated[i, j] += img[i - y_tc, j - x_tc] * x_p * y_p
                translated[i, j] += img[i - y_tf, j - x_tc] * x_p * (1 - y_p)
                translated[i, j] += img[i - y_tc, j - x_tf] * (1 - x_p) * y_p
                translated[i, j] += img[i - y_tf, j - x_tf] * (1 - x_p) * (1 - y_p)

    return translated

ch = 0

hybs = 80
#stack = np.zeros([2048,2048,hybs+2])
stack = np.zeros([hybs+2,2048,2048])


#stack[:,:,0] = skimage.io.imread(snakemake.input[1])
stack[0,:,:] = skimage.io.imread(snakemake.input[1])

offsets = pd.read_csv(snakemake.input[0])
offsets.loc[:,"x"] = -offsets.loc[:,"x"]
offsets.loc[:,"y"] = -offsets.loc[:,"y"]

for i in range(2, len(snakemake.input)):
    print("hyb ", i-1)
    img = skimage.io.imread(snakemake.input[i])
    trans_img = translate_img(img, offsets, i-2)
    #stack[:,:,i-1] = trans_img
    stack[i-1,:,:] = trans_img

tifffile.imsave(snakemake.output[0], stack)
