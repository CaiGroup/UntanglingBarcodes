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

"""
ch = 0
pos = 0

min_params = 'lf80.0_wf12.0_sf0.0_dr0'

pnts_fname = "../results/points_ch_{ch}_pos_{pos}_".format(ch=ch, pos=pos)+min_params+"_mw600_rbr_3.3.csv"
#aligned_stack_fname = "../alignment/aligned_stack_ch_{ch}_pos_{pos}.tif".format(ch=ch, pos=pos)
aligned_stack_fname = "../plots/aligned_stacks/aligned_stack_ch_{ch}_pos_{pos}.tif".format(ch=ch, pos=pos)
published_pnts_fname = "../validation_files/published_pnt_locs_run2.csv"
#gene_cb_fname = "../validation_files/10k_codebook_ch_{ch}.csv".format(ch=ch)
gene_cb_fname = "../codebooks/10k_codebook_ch_{ch}.csv".format(ch=ch)
full_cb_fname = "../codebooks/codebook_ch_{ch}.csv".format(ch=ch)
roi_width = 200

"""

ch = int(snakemake.wildcards['ch'])
pos = int(snakemake.wildcards['pos'])

pnts_fname = snakemake.input[0]
aligned_stack_fname = snakemake.input[1]
published_pnts_fname = snakemake.input[2]
gene_cb_fname = snakemake.input[3]
full_cb_fname = snakemake.input[4]
#rois = snakemake.imput[5]
roi_width = snakemake.params["roi_width"]


aligned_stack_ = skimage.io.imread(aligned_stack_fname)#"aligned_stack_ch_2.tif")


plot_width = 800
plot_height = 800
#xstart = 1000
#ystart = 1000
ystart, xstart = auto_roi(aligned_stack_, roi_width)
xend = xstart + roi_width#1200
yend = ystart + roi_width#1200

pnts = pd.read_csv(pnts_fname)#'results/points_ch_2_lf20.0_wf8.0_sf0.0_dr0_mw1000_rbr_3.3.csv')
pnts.x = pnts.x -2
pnts.y = pnts.y -2
pnts['hybridization'] = pnts['hyb']


#aligned_stack = aligned_stack_[1:81, xstart:xend, ystart:yend]
#aligned_stack = aligned_stack_[0:80, xstart:xend, ystart:yend]
aligned_stack = aligned_stack_[1:81, ystart:yend, xstart:xend]

n_ims, height, width = np.shape(aligned_stack)
ashv = hv.Dataset((np.arange(width), np.arange(height), np.arange(1, 81), aligned_stack),
                ['x', 'y', 'hybridization'], 'Fluorescence')

def get_roi_inds(df):
    return np.logical_and(np.logical_and(df.x > xstart, df.y > ystart), np.logical_and(df.x < xend, df.y < yend))

pnts = pnts.loc[get_roi_inds(pnts), :]#["hybridization", "x", "y", "result"]]
pnts['x'] = pnts['x'] - xstart
pnts['y'] = pnts['y'] - ystart

publ_points = pd.read_csv(published_pnts_fname)#"published_pnt_locs_run2.csv")
pub_points = publ_points.loc[publ_points['pos'] == pos]
pub_points = pub_points.loc[get_roi_inds(pub_points)]
#cb = pd.read_csv("validation_files/10k_barcodes_488.csv",names = ["gene_name", "fpkm?", '1', '2', '3', '4'])
cb = pd.read_csv(gene_cb_fname)#,names = ["gene_name", '1', '2', '3', '4'])

#cb = pd.read_csv("validation_files/10k_647_ALL_3334.csv",names = ["gene_name", '1', '2', '3', '4'])
pub_points_w_hybs = pub_points.merge(cb, on="gene_name")
pub_points_w_hybs['codeword'] = ['%s, %s, %s, %s' % tuple(pub_points_w_hybs.loc[i,['Round1','Round2','Round3','Round4']]) for i in range(len(pub_points_w_hybs))]
#pub_points_w_hybs.drop("Unnamed: 0", axis='columns', inplace=True)
#pub_points_w_hybs.drop('fpkm?', axis='columns', inplace=True)
id_vars = ['pos', 'cellid', 'gene_name', 'gene_num', 'x', 'y', 'codeword']
pub_pointsmt = pub_points_w_hybs.melt(id_vars=id_vars, var_name='round',value_name="pseudocolor")
pub_pointsmt["round"] = [int(rnd[-1]) for rnd in pub_pointsmt["round"]]
pub_pointsmt['pseudocolor'] = [int(pc) for pc in pub_pointsmt["pseudocolor"]]
pub_pointsmt['hybridization'] = 20*(pub_pointsmt['round']-1) + pub_pointsmt['pseudocolor'] - 1
pps = pub_pointsmt.loc[:,['hybridization', 'x', 'y', 'gene_name', 'codeword']]
pps['x'] = pps['x'] - xstart -1
pps['y'] = pps['y'] - ystart -1
pps['hybridization'] = pps['hybridization']+1

pnts['gene_number'] = pnts['decoded']
cb["gene_number"] = np.array([int(ind) for ind in cb.index]) + 1
on_targets = pnts.merge(cb, on='gene_number')
on_targets['codeword'] = ['%s, %s, %s, %s' % tuple(on_targets.loc[i,['Round1','Round2','Round3','Round4']]) for i in range(len(on_targets))]


#full_cb = pd.read_table('codebooks/E2019_cb_all_control_488.txt', names=['r1','r2','r3','r4'])
full_cb = pd.read_csv(full_cb_fname)#, names=['r1','r2','r3','r4'])

full_cb['gene_number'] = range(1,len(full_cb)+1)
off_targets = pnts[pnts.decoded > max(cb.gene_number)]
off_targets = off_targets.merge(full_cb, on='gene_number')
off_targets.reset_index(inplace=True)
off_targets['codeword'] = ['%s, %s, %s, %s' % tuple(off_targets.loc[i,['Round1','Round2','Round3','Round4']]) for i in range(len(off_targets))]

published_tool_tips = [
    ('Gene Name', '@gene_name'),
    ('Codeword', '@codeword')
]

ontarget_tool_tips = [
    ('Gene Name', '@gene_name'),
    ('Codeword', '@codeword'),
    ('Weight', '@w'),
    ('Sigma', '@s')
]

offtarget_tool_tips = [
    ('Codeword', '@codeword'),
    ('Weight', '@w'),
    ('Sigma', '@s')
]

nd_tool_tips = [
    ('Weight', '@w'),
    ('Sigma', '@s')
]

published_hover = HoverTool(tooltips=published_tool_tips)
ontarget_hover = HoverTool(tooltips=ontarget_tool_tips)
offtarget_hover = HoverTool(tooltips=offtarget_tool_tips)
nd_hover = HoverTool(tooltips=nd_tool_tips)

ppsds = hv.Dataset(pps, ['hybridization', 'x', 'y'], ['gene', 'codeword'])
ontgt = hv.Dataset(on_targets, ['hybridization', 'x', 'y'], ['gene', 'codeword', 's', 'w'])
offtgt = hv.Dataset(off_targets, ['hybridization', 'x', 'y'], ['codeword', 's', 'w'])
nd = hv.Dataset(pnts[pnts.decoded == 0], ['hybridization', 'x', 'y'], ['s', 'w'])
ims = ashv.to(hv.Image, ['x', 'y'])


pps_hv = ppsds.to(hv.Points, ['x', 'y'], label='Published')
ontarget_hv = ontgt.to(hv.Points, ['x', 'y'], label='On-Target')
offtarget_hv = offtgt.to(hv.Points, ['x', 'y'], label='Off-Target')
nd_hv = nd.to(hv.Points, ['x', 'y'], label='Not-Decoded')

pps_hv.opts(tools=[published_hover], size=5,color='black', marker='*', width = plot_width, height=plot_height)
ontarget_hv.opts(tools=[ontarget_hover], size=4,color='blue', width = plot_width, height=plot_height)
offtarget_hv.opts(tools=[offtarget_hover], size=5,color='red', marker='x', width = plot_width, height=plot_height)
nd_hv.opts(tools=[nd_hover], size=5, color='green', marker='+', width = plot_width, height=plot_height)
ims.opts(cmap='viridis', colorbar=True, width=plot_width, height=plot_height)

overlay = ims * pps_hv * ontarget_hv * offtarget_hv * nd_hv

hv.save(overlay, snakemake.output[0])
#hv.save(overlay, 'test_overlay.html')
