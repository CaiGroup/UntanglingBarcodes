import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np

# find barcodes/gene for each cell in fov
ninputs = len(snakemake.input)

print("input files:")
for file in snakemake.input:
	print(file)

def count_cells(df):
    return len(list(df.groupby(['pos','cellid'])))

results_2019_rep_array = [pd.read_csv(snakemake.input[i]) for i in [0,1]]
ncells = sum(map(count_cells, results_2019_rep_array))
results_2019 = pd.concat(results_2019_rep_array)

results_2022 = pd.concat([pd.read_csv(snakemake.input[i]) for i in [2,3]])
#results_2019 = pd.read_csv("resources/validation_files/published_pnt_locs_run2.csv")
#results_new = pd.read_csv("results/cell_barcode_counts/cell_barcode_counts_lf30.0_wf4.0_sf0.0_dr0.csv")
smFISH_fname = snakemake.input[4]


smFISH_results = pd.read_csv(smFISH_fname)
ncells_smFISH = len(smFISH_results)
smFISH_counts = np.sum(smFISH_results, axis=0)/ncells_smFISH
print("smFISH head")
print(smFISH_counts.head())


def process_2022_results(results_2022):

    results_2022 = results_2022.loc[results_2022.gene != "negative_control"]
    results_2022_g = results_2022.groupby(['ch', 'gene']).sum()#.loc["counts"]
    results_2022.name= "new_counts"
    #results_new.set_index('gene', inplace=True)

    results_2022_g =results_2022_g.reset_index()
    results_2022_g.set_index('gene', inplace=True)
    return results_2022_g


#results_2019.rename(columns={'fov':'pos', 'cell':'cellid'}, inplace=True)

results_2019_cnts = results_2019.groupby("gene").size()

results_new_g = process_2022_results(results_2022)

print("results_new_g.head()")
print(results_new_g.head())

print("ncells: ", ncells)

results_new_g['old_counts'] = results_2019_cnts

results_new_g = results_new_g.loc[:,["ch","counts","old_counts"]]
results_new_g= results_new_g.reset_index()
channels = np.unique(results_new_g.ch)
results_new_g.set_index(["ch","gene"], inplace=True)

results_new_g.rename(columns={'counts':'new_counts'}, inplace=True)

results_new_g /= ncells

for i, ch in enumerate(channels):

    lm_synd_pub = stats.linregress(results_new_g.loc[ch, 'old_counts'], results_new_g.loc[ch,'new_counts'])
    title = 'channel %d: r = %.2f, m = %.2f, b = %.2f' % (ch, lm_synd_pub.rvalue, lm_synd_pub.slope, lm_synd_pub.intercept)

    results_new_g.loc[ch].plot("old_counts", 'new_counts','scatter',title=title)
    plt.savefig(snakemake.output[i])

#results_new_g.loc[488].plot("old_counts", 'counts','scatter',title="488")
#results_new_g.loc[561].plot("old_counts", 'counts','scatter',title="561")
#results_new_g.loc[643].plot("old_counts", 'counts','scatter',title="643")

#df = pd.DataFrame({"old":results_2019_cnts, "new": results_new_cnts["counts"], "ch":results_new_cnts["ch"]})
#results_2019_cnts.loc['old'] = results_2019_cnts

#results_2019_cnts.fillna(0, inplace=True)


#df /=ncells

#df_488 = df.loc[df]


lm_synd_pub = stats.linregress(results_new_g['old_counts'], results_new_g['new_counts'])

title = 'All Channels: r = %.2f, m = %.2f, b = %.2f' % (lm_synd_pub.rvalue, lm_synd_pub.slope, lm_synd_pub.intercept)

results_new_g.plot("old_counts", "new_counts", "scatter", title=title)
plt.savefig(snakemake.output[3])
plt.clf()

print("results_new_g.index")
print(results_new_g.index)
results_new_g.reset_index("ch", drop=True, inplace=True)

smFISH_counts.index.name='gene'
print("smFISH head")
print(smFISH_counts.head())
results_new_g['smFISH'] = smFISH_counts

print("results with smFISH")
print(results_new_g.head())

results_new_g.dropna(inplace=True)
print("results with smFISH, dropna")
print(results_new_g.head())

smFISH_new_fit = stats.linregress(results_new_g['smFISH'], results_new_g['new_counts'])
smFISH_pub_fit = stats.linregress(results_new_g['smFISH'], results_new_g['old_counts'])

plt.plot(results_new_g['smFISH'], results_new_g['old_counts'], '.', label="2019; m=%.2f, r=%.2f" % (smFISH_pub_fit.slope, smFISH_pub_fit.rvalue))
plt.plot(results_new_g['smFISH'], results_new_g['new_counts'], '.', label="2022; m=%.2f, r=%.2f" % (smFISH_new_fit.slope, smFISH_new_fit.rvalue))
plt.xlabel("smFISH counts")
plt.ylabel("seqFISH+ counts")
plt.legend()
plt.savefig(snakemake.output[4])


"""
rnaseq = pd.read_csv("resources/validation_files/nih3t3_FPKM.csv")
rnaseq.rename(columns={'tracking_id':'gene'}, inplace=True)
rnaseq.dropna(inplace=True)
rnaseq.set_index("gene",inplace=True)
stats.linregress(rnaseq['new'], rnaseq['3T3 B1'])
rnaseq.plot('new', '3T3 B1', 'scatter',logx=True,logy=True,s=0.1)
plt.savefig(snakemake.output[4])
"""

"""
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
"""
