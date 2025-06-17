# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 14:30:07 2021

@author: jonat
"""
import skimage.io
import holoviews as hv
hv.extension('bokeh')
import numpy as np
import pandas as pd
from bokeh.models import HoverTool
from auto_roi import auto_roi

ch = int(snakemake.wildcards['ch'])
fov = int(snakemake.wildcards['pos'])
pnts_fname = snakemake.input[0]
im_stack_fname = snakemake.input[1]
roi_width = snakemake.params["roi_width"]

im_stack_ = skimage.io.imread(im_stack_fname)#"aligned_stack_ch_2.tif")

plot_width = 800
plot_height = 800

ystart, xstart = auto_roi(im_stack_, roi_width)
xend = xstart + roi_width
yend = ystart + roi_width

pnts = pd.read_csv(pnts_fname)
pnts.x = pnts.x -1
pnts.y = pnts.y -1
pnts.z = pnts.z -1

im_stack = im_stack_[:, ystart:yend, xstart:xend]

depth, height, width = np.shape(im_stack)
ashv = hv.Dataset((np.arange(width), np.arange(height), np.arange(depth), im_stack),
                ['x', 'y', 'z'], 'Fluorescence')

def get_roi_inds(df):
    return np.logical_and(np.logical_and(df.x > xstart, df.y > ystart), np.logical_and(df.x < xend, df.y < yend))

pnts = pnts.loc[get_roi_inds(pnts), :]#["hybridization", "x", "y", "result"]]
pnts.loc[:,"z"] = np.round(pnts.z)
pnts.loc[pnts.z < 0,"z"] = 0
pnts.loc[:,"y"] -= ystart
pnts.loc[:,"x"] -= xstart

tool_tips = [
    ('Weight', '@w'),
    ('Sigma_xy', '@s_xy'),
    ('Sigma_z', '@s_z')
]

hover = HoverTool(tooltips=tool_tips)

pnts_ds = hv.Dataset(pnts, ['x', 'y', 'z'], ['sxy', 'sz', 'w'])
ims = ashv.to(hv.Image, ['x', 'y'])

pnts_hv = pnts_ds.to(hv.Points, ['x', 'y'])

pnts_hv.opts(tools=[hover], size=4,color='blue', width = plot_width, height=plot_height)
ims.opts(cmap='viridis', colorbar=True, width=plot_width, height=plot_height)

overlay = ims * pnts_hv

hv.save(overlay, snakemake.output[0])
