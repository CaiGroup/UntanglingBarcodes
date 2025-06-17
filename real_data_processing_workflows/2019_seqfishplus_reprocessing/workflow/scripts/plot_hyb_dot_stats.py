import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

cb_fname = snakemake.input[0]
cb=pd.read_csv(cb_fname)

pnts_fname = snakemake.input[1]
ch = pnts_fname.split('_')[2]
pos = pnts_fname.split('_')[4]
ps = pd.read_csv(pnts_fname)
nd = ps.loc[ps.decoded==0]
max_on_target_cw = sum(cb['gene_name'] != 'negative_control')
off_target = ps.loc[ps.decoded > max_on_target_cw]
on_target = ps.loc[np.logical_and(ps.decoded <= max_on_target_cw, ps.decoded != 0)]


ndg = nd.groupby('hyb')
on_targetg = on_target.groupby('hyb')
oftg = off_target.groupby('hyb')
ag = ps.groupby('hyb')

ws = pd.DataFrame({'all':ag.mean()['w'],
                   'off target':oftg.mean()['w'],
                   'on target':on_targetg.mean()['w'],
                   'not decoded':ndg.mean()['w']})
ws.plot()
plt.xlabel("Hybridization")
plt.ylabel("Average Dot Intensity (a.u.)")
plt.title('Position {pos} channel {ch}'.format(pos=pos,ch=ch))
plt.savefig(snakemake.output[0])

szs = pd.DataFrame({'all':ag.size(),
                   'off target':oftg.size(),
                   'on target':on_targetg.size(),
                   'not decoded':ndg.size()})
szs.plot()
plt.xlabel("Hybridization")
plt.ylabel("Number of Dots")
plt.title('Position {pos} channel {ch}'.format(pos=pos,ch=ch))
plt.savefig(snakemake.output[1])


"""
sns.ecdfplot(on_target.set_index('hyb').loc[1], x='w', label="Decoded")
sns.ecdfplot(nd.set_index('hyb').loc[1], x='w', label="Not Decoded")
sns.ecdfplot(off_target.set_index('hyb').loc[1], x='w', label="Off Target")


plt.legend()
plt.xlabel("Intensity (a.u.)")
plt.title('ECDF Position {pos} channel {ch}'.format(pos=pos,ch=ch))
plt.savefig(snakemake.output[2])
"""
