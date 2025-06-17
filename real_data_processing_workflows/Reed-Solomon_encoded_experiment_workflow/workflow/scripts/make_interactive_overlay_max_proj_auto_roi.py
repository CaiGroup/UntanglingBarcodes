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
from scipy.ndimage import shift

ch = int(snakemake.wildcards['ch'])
fov = int(snakemake.wildcards['pos'])
pnts_fname = snakemake.input[0]
offsets_fname = snakemake.input[1]

roi_width = snakemake.params["roi_width"]

offsets = pd.read_csv(offsets_fname)

hyb_stack = np.zeros([100,2048,2048])
for hyb, im_stack_fname in enumerate(snakemake.input[2:]):
#for hyb in range(36):
    im_stack = skimage.io.imread(im_stack_fname)#"aligned_stack_ch_2.tif")
    
    #im_stack = skimage.io.imread('../results/ch_stacks/HybCycle_%d_ch_%d_pos%d.tif' % (hyb, ch, fov))
    max_proj = im_stack.max(axis=0)
    hyb_stack[hyb, :, :] = shift(max_proj, (-offsets.loc[hyb,'y_offset'], -offsets.loc[hyb,'x_offset']), order=1)
    
    


plot_width = 800
plot_height = 800

ystart, xstart = auto_roi(hyb_stack, roi_width)
xend = xstart + roi_width
yend = ystart + roi_width

pnts = pd.read_csv(pnts_fname)
pnts.x = pnts.x - 1
pnts.y = pnts.y - 1
pnts.hyb = pnts.hyb - 1

hyb_stack = hyb_stack[:, ystart:yend, xstart:xend]

depth, height, width = np.shape(hyb_stack)
ashv = hv.Dataset((np.arange(width), np.arange(height), np.arange(depth), hyb_stack),
                ['x', 'y', 'hyb'], 'Fluorescence')

def get_roi_inds(df):
    return np.logical_and(np.logical_and(df.x > xstart, df.y > ystart), np.logical_and(df.x < xend, df.y < yend))

pnts = pnts.loc[get_roi_inds(pnts), :]#["hybridization", "x", "y", "result"]]
pnts.loc[:,"z"] = np.round(pnts.z)
pnts.loc[:,"y"] -= ystart
pnts.loc[:,"x"] -= xstart

tool_tips = [
    ('Weight', '@w'),
    ('Sigma_xy', '@s_xy'),
    ('Sigma_z', '@s_z')
]

hover = HoverTool(tooltips=tool_tips)

undecoded_pnts = pnts.loc[pnts.decoded == 0, :]
offtarget_pnts = pnts.loc[pnts.decoded > 581]
ontarget_pnts = pnts.loc[(pnts.decoded < 582) & (pnts.decoded != 0)]



undecoded_pnts_ds = hv.Dataset(undecoded_pnts, ['x', 'y', 'hyb'], ['sxy', 'sz', 'w'])
ontarget_pnts_ds = hv.Dataset(ontarget_pnts, ['x', 'y', 'hyb'], ['sxy', 'sz', 'w'])
offtarget_pnts_ds = hv.Dataset(offtarget_pnts, ['x', 'y', 'hyb'], ['sxy', 'sz', 'w'])
#ims = ashv.to(hv.Image, ['x', 'y'])

undecoded_pnts_hv = undecoded_pnts_ds.to(hv.Points, ['x', 'y'])
ontarget_pnts_hv = ontarget_pnts_ds.to(hv.Points, ['x', 'y'])
offtarget_pnts_hv = offtarget_pnts_ds.to(hv.Points, ['x', 'y'])

undecoded_pnts_hv.opts(tools=[hover], size=4,color='blue', width = plot_width, height=plot_height)
ontarget_pnts_hv.opts(tools=[hover], size=4,color='black', width = plot_width, height=plot_height)
offtarget_pnts_hv.opts(tools=[hover], size=4,color='red', width = plot_width, height=plot_height)
#ims.opts(cmap='viridis', colorbar=True, width=plot_width, height=plot_height)


ims = ashv.to(hv.Image, ['x', 'y'])


ims.opts(cmap='viridis', colorbar=True, width=plot_width, height=plot_height)

overlay = ims * undecoded_pnts_hv * ontarget_pnts_hv * offtarget_pnts_hv

hv.save(overlay, snakemake.output[0])
