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
#from bokeh.models import HoverTool



#ch = int(snakemake.wildcards['ch'])
#pos = int(snakemake.wildcards['pos'])

true_pnts_fname = snakemake.input[0]
processed_pnts_fname = snakemake.input[1]
full_cb_fname = snakemake.input[2]
#rois = snakemake.imput[5]
fov_width = snakemake.params["fov_width"]



#stack = np.zeros([fov_width,fov_width, 80])
stack = np.zeros([80, fov_width,fov_width])


for i in range(3, len(snakemake.input)):
    print("hyb ", i-2)
    #stack[:,:,i-3] = skimage.io.imread(snakemake.input[i])
    stack[i-3,:,:] = skimage.io.imread(snakemake.input[i])



true_pnts = pd.read_csv(true_pnts_fname)#'results/points_ch_2_lf20.0_wf8.0_sf0.0_dr0_mw1000_rbr_3.3.csv')
true_pnts.x = true_pnts.x -2
true_pnts.y = true_pnts.y -2
true_pnts['hybridization'] = true_pnts['hyb']


#aligned_stack = aligned_stack_[1:81, xstart:xend, ystart:yend]
#aligned_stack = aligned_stack_[0:80, xstart:xend, ystart:yend]

n_ims = 80
ashv = hv.Dataset((np.arange(fov_width), np.arange(fov_width), np.arange(1, 81), stack),
                ['x', 'y', 'hybridization'], 'Fluorescence')

"""
pnts = pnts.loc[get_roi_inds(pnts), :]#["hybridization", "x", "y", "result"]]
pnts['x'] = pnts['x'] - xstart
pnts['y'] = pnts['y'] - ystart

publ_points = pd.read_csv(published_pnts_fname)#"published_pnt_locs_run2.csv")
#cb = pd.read_csv("validation_files/10k_barcodes_488.csv",names = ["gene_name", "fpkm?", '1', '2', '3', '4'])
cb = pd.read_csv(gene_cb_fname)#,names = ["gene_name", '1', '2', '3', '4'])

#cb = pd.read_csv("validation_files/10k_647_ALL_3334.csv",names = ["gene_name", '1', '2', '3', '4'])
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

ppsds = hv.Dataset(pps, ['hybridization', 'x', 'y'], ['gene_name', 'codeword'])
ontgt = hv.Dataset(on_targets, ['hybridization', 'x', 'y'], ['gene_name', 'codeword', 's', 'w'])
offtgt = hv.Dataset(off_targets, ['hybridization', 'x', 'y'], ['codeword', 's', 'w'])
nd = hv.Dataset(pnts[pnts.decoded == 0], ['hybridization', 'x', 'y'], ['s', 'w'])
"""
true_pnts_tool_tips = [
    ('Gene Number', '@gene_number')
]

ims = ashv.to(hv.Image, ['x', 'y'])
print(true_pnts.head())
true_pnts_ds  = hv.Dataset(true_pnts, ['hybridization', 'x', 'y'], ['gene_number'])

true_pnts_hv = true_pnts_ds.to(hv.Points, ['x', 'y'], label='True')

#pps_hv = ppsds.to(hv.Points, ['x', 'y'], label='Published')
#ontarget_hv = ontgt.to(hv.Points, ['x', 'y'], label='On-Target')
#offtarget_hv = offtgt.to(hv.Points, ['x', 'y'], label='Off-Target')
#nd_hv = nd.to(hv.Points, ['x', 'y'], label='Not-Decoded')

#true_pnts_hv.opts(tools=[true_pnts_tool_tips], size=5,color='black', marker='*')#, width = fov_width, height=fov_width)

true_pnts_hv.opts(size=5,color='black', marker='*') #, width = fov_width, height=fov_width)

#ontarget_hv.opts(tools=[ontarget_hover], size=4,color='blue', width = plot_width, height=plot_height)
#offtarget_hv.opts(tools=[offtarget_hover], size=5,color='red', marker='x', width = plot_width, height=plot_height)
#nd_hv.opts(tools=[nd_hover], size=5, color='green', marker='+', width = plot_width, height=plot_height)
ims.opts(cmap='viridis')#, colorbar=True) #, width=fov_width, height=fov_width)

overlay = ims  * true_pnts_hv

hv.save(overlay, snakemake.output[0])
#hv.save(overlay, 'test_overlay.html')
