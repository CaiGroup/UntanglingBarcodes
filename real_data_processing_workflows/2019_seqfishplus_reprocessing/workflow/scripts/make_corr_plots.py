# -*- coding: utf-8 -*-
"""
Created on Thu Aug 12 15:36:14 2021

@author: jonat
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from get_segmented_cells import get_segmented_cells
import scipy.stats as stats
import skimage.io

mpaths_fname = snakemake.input[0]
published_pnts_fname = snakemake.input[1]
smFISH_fname = snakemake.input[2]
gene_cb_fname = snakemake.input[3]
labeled_img_name = snakemake.input[4]
pnts_fname = snakemake.input[5]

#ch_wvln = {'0':647,'1':561,'2':488}
wvln = int(snakemake.wildcards["ch"])

mpaths_df = pd.read_csv(mpaths_fname)#"results/mpaths_ch_2_lf20.0_wf8.0_sf0.0_dr0_mw1000.csv")
publ_points = pd.read_csv(published_pnts_fname)#"../published_pnt_locs_run2.csv")
publ_points = publ_points.loc[publ_points['pos'] == int(snakemake.wildcards['pos'])]

#mpaths_df = pd.read_csv("20210528/20210528_results/mpaths_lf40.0_wf12.0_sf0.0_mw600.csv")
#mpaths_df = pd.read_csv("results/mpaths_ch_1_lf40.0_wf8.0_sf0.0_dr0_mw600.csv")

aligned_points = pd.read_csv(pnts_fname)

"""
if snakemake.wildcards["ch"] == '0':
    ontarget_cutoff = 3335
else:
    ontarget_cutoff = 3334
"""
# keep only on-target barcodes
#mpaths_df = mpaths_df.loc[mpaths_df.gene_number < ontarget_cutoff]
mpaths_df = mpaths_df.loc[mpaths_df.gene != "negative_control"]

cb = pd.read_csv(gene_cb_fname)

smFISH_results = pd.read_csv(smFISH_fname)
ncells_smFISH = len(smFISH_results)
smFISH_counts = np.sum(smFISH_results, axis=0)/ncells_smFISH

sd_gene_cnts = mpaths_df.groupby('gene').size()

pub_gene_cnts = publ_points.groupby('gene').size()



res  = cb.drop(['block1', 'block2', 'block3', 'block4'], axis=1)
res.set_index('gene', inplace=True)
res['new'] = sd_gene_cnts
res['published'] = pub_gene_cnts

res.fillna(0, inplace = True)

labeled_img = skimage.io.imread(labeled_img_name)
#regions = np.unique(labeled_img)
#print("regions: ", regions)
n_cells = len(list(filter(lambda x : x != 0, np.unique(labeled_img))))
#print("n cells: ", n_cells)
res = res/n_cells
#res = res/len(rois)



smFISH_counts.index.name='gene'
res['smFISH'] = smFISH_counts


lm_synd_pub = stats.linregress(res['published'], res['new'])

r2 = res.dropna()
smFISH_new_fit = stats.linregress(r2['smFISH'], r2['new'])

smFISH_pub_fit = stats.linregress(r2['smFISH'], r2['published'])


fig = plt.figure(figsize = (10,3))
#fig, ax = plt.subplots(1,3)
ax = fig.subplots(1,3, sharey=True)



ax[0].plot(res['published'], res['new'], '.')
ax[0].set_xlabel('Gene Trancript Counts Published')
ax[0].set_ylabel('Gene Transcript Counts New')
ax[0].set_title('r = %.2f, m = %.2f, b = %.2f' % (lm_synd_pub.rvalue, lm_synd_pub.slope, lm_synd_pub.intercept))

ax[1].plot(r2['smFISH'], r2['published'], '.')

ax[1].set_xlabel('Gene Transcript Counts smFISH')
ax[1].set_ylabel('Gene Trancript Counts Published')
ax[1].set_title('r = %.2f, m = %.2f, b = %.2f' % (smFISH_pub_fit.rvalue, smFISH_pub_fit.slope, smFISH_pub_fit.intercept))



ax[2].plot(r2['smFISH'], r2['new'], '.')

ax[2].set_xlabel('Gene Transcript Counts smFISH')
ax[2].set_ylabel('Gene Transcript Counts New')
#ax[2].title('Published vs New Counts Channel ')
ax[2].set_title('r = %.2f, m = %.2f, b = %.2f' % (smFISH_new_fit.rvalue, smFISH_new_fit.slope, smFISH_new_fit.intercept))

plt.tight_layout()
plt.savefig(snakemake.output[0])

res.to_csv(snakemake.output[1])
r2.to_csv(snakemake.output[2])


### Filtered
"""
mpaths_df = pd.read_csv(mpaths_fname)

def neg_cntrl_saturated_ccs(df):
    n_nctrls = np.sum(df.gene == "negative_control")
    remove_if_true = len(df) > 5 and n_nctrls/len(df) > 0.3
    #print("remove_if_true: ",  remove_if_true)
    return not remove_if_true

mpaths_df_filtered = mpaths_df.groupby("cc").filter(neg_cntrl_saturated_ccs)

sd_gene_cnts = mpaths_df_filtered.groupby('gene').size()

pub_gene_cnts = publ_points.groupby('gene').size()



res  = cb.drop(['block1', 'block2', 'block3', 'block4'], axis=1)
res.set_index('gene', inplace=True)
res['new'] = sd_gene_cnts
res['published'] = pub_gene_cnts

res.fillna(0, inplace = True)

labeled_img = skimage.io.imread(labeled_img_name)
#regions = np.unique(labeled_img)
#print("regions: ", regions)
n_cells = len(list(filter(lambda x : x != 0, np.unique(labeled_img))))
#print("n cells: ", n_cells)
res = res/n_cells
#res = res/len(rois)



smFISH_counts.index.name='gene'
res['smFISH'] = smFISH_counts


lm_synd_pub = stats.linregress(res['published'], res['new'])

r2 = res.dropna()
smFISH_new_fit = stats.linregress(r2['smFISH'], r2['new'])

smFISH_pub_fit = stats.linregress(r2['smFISH'], r2['published'])


fig = plt.figure(figsize = (10,3))
#fig, ax = plt.subplots(1,3)
ax = fig.subplots(1,3, sharey=True)



ax[0].plot(res['published'], res['new'], '.')
ax[0].set_xlabel('Gene Trancript Counts Published')
ax[0].set_ylabel('Gene Transcript Counts New')
ax[0].set_title('r = %.2f, m = %.2f, b = %.2f' % (lm_synd_pub.rvalue, lm_synd_pub.slope, lm_synd_pub.intercept))

ax[1].plot(r2['smFISH'], r2['published'], '.')

ax[1].set_xlabel('Gene Transcript Counts smFISH')
ax[1].set_ylabel('Gene Trancript Counts Published')
ax[1].set_title('r = %.2f, m = %.2f, b = %.2f' % (smFISH_pub_fit.rvalue, smFISH_pub_fit.slope, smFISH_pub_fit.intercept))



ax[2].plot(r2['smFISH'], r2['new'], '.')

ax[2].set_xlabel('Gene Transcript Counts smFISH')
ax[2].set_ylabel('Gene Transcript Counts New')
#ax[2].title('Published vs New Counts Channel ')
ax[2].set_title('r = %.2f, m = %.2f, b = %.2f' % (smFISH_new_fit.rvalue, smFISH_new_fit.slope, smFISH_new_fit.intercept))

plt.tight_layout()
plt.savefig(snakemake.output[3])
"""
