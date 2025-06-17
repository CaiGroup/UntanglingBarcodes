# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 11:48:25 2021

@author: jonat
"""
from util import pil_imread
import xarray
import numpy as np
import tifffile
import skimage
import skimage.io

#im_stack = skimage.io.imread(snakemake.input[0])
def get_ch_im(snakemake, tifstack_name):
    if "z" in snakemake.wildcards.keys():
        z = int(snakemake.wildcards["z"])
    else:
        z = snakemake.params["z"]

    if "ch" in snakemake.wildcards.keys():
        save_ch = int(snakemake.wildcards["ch"])
    else:
        save_ch = snakemake.params["ch"]


    tifstack, metadata = pil_imread(tifstack_name,metadata=True)

    channels = np.unique([int(immddict['Channel']) for immddict in metadata])
    channels=-np.sort(-channels)
    nchannels, nzs, nys, nxs = np.shape(tifstack)
    ds = xarray.Dataset({'fluorescence': (("channel", "z", "y", "x"), tifstack)},
                             {"channel":channels,
                              "z": range(nzs),
                              "y": range(nys),
                              "x": range(nxs)
                             })

    ch_im = ds.sel(channel=save_ch, z=z)
    return ch_im


for i, infile in enumerate(snakemake.input):
    #im_stack = skimage.io.imread(infile)
    hyb = int(infile.split('_')[1].split('/')[0])
    if hyb == 61:
        assert(int(snakemake.output[i+1].split('_')[2]) == 62)
        ch_im = get_ch_im(snakemake, infile)
        #tifffile.imsave(snakemake.output[i+1], ch_im.fluorescence)
        skimage.io.imsave(snakemake.output[i+1], ch_im.fluorescence)
    elif hyb == 62:
        assert(int(snakemake.output[i-1].split('_')[2]) == 61)
        ch_im = get_ch_im(snakemake, infile)
        #tifffile.imsave(snakemake.output[i-1], ch_im.fluorescence)
        skimage.io.imsave(snakemake.output[i-1], ch_im.fluorescence)
    else:
        assert(int(snakemake.output[i].split('_')[2]) == hyb)
        ch_im = get_ch_im(snakemake, infile)
        #tifffile.imsave(snakemake.output[i], ch_im.fluorescence)
        skimage.io.imsave(snakemake.output[i], ch_im.fluorescence)
