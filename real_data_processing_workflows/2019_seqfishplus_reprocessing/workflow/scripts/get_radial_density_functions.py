import pandas as pd
import numpy as np
from scipy.spatial import cKDTree as KDTree
import matplotlib.pyplot as plt

nchannels=snakemake.params['nchannels']
r_max = snakemake.params['rmax']
delta_r = snakemake.params['delta_r']

def weighted_radial_density(points, rs):
    ps_tree = KDTree(points.loc[:,['x','y']])
    ps_nbrs = ps_tree.count_neighbors(ps_tree, rs)
    nself_ps = ps_nbrs - len(points)
    return (nself_ps[1:] - nself_ps[:-1])/(np.pi*(rs[1:]**2 -rs[:-1]**2))


#dt_cnts = pd.DataFrame(columns=['all','on_target', 'off_target', 'not_decoded'])
dt_cnts = pd.DataFrame(np.zeros([1,4], np.int64), columns=['all','on_target', 'off_target', 'not_decoded'])

pnts_fnames = snakemake.input[nchannels:]
cb_fnames = snakemake.input[:nchannels]
max_on_target_cw = {}
for cb_fname in cb_fnames:
    ch = ''.join(filter(lambda c: c.isdigit(), cb_fname))
    cb=pd.read_csv(cb_fname)
    max_on_target_cw[ch] = sum(cb['gene'] != 'negative_control')

rs = np.arange(0, r_max, delta_r)

weighted_densities = pd.DataFrame(np.zeros([len(rs)-1,4]),columns=['all','on_target','off_target','not_decoded'])

for pnts_fname in snakemake.input[nchannels:]:
    print(pnts_fname)
    #ch = pnts_fname.split('_')[2]
    ch = pnts_fname.split('_')[5]
    print("ch: ", ch)
    ps = pd.read_csv(pnts_fname)
    nd = ps.loc[ps.decoded==0]
    off_target = ps.loc[ps.decoded > max_on_target_cw[ch]]
    on_target = ps.loc[np.logical_and(ps.decoded <= max_on_target_cw[ch], ps.decoded != 0)]

    weighted_densities.loc[:,'all'] += weighted_radial_density(ps, rs)
    weighted_densities.loc[:,'not_decoded'] += weighted_radial_density(nd, rs)
    weighted_densities.loc[:,'off_target'] += weighted_radial_density(off_target, rs)
    weighted_densities.loc[:,'on_target'] += weighted_radial_density(on_target, rs)

    dt_cnts+=np.array([len(ps), len(on_target), len(off_target), len(nd)])

print(dt_cnts)
print(weighted_densities.head())
radial_densities = np.divide(weighted_densities,dt_cnts)
radial_densities.to_csv(snakemake.output[0], index=False)
print(radial_densities.head())
plt.plot(rs[1:], radial_densities['all'], label="all")
plt.plot(rs[1:], radial_densities['not_decoded'], label="not decoded")
plt.plot(rs[1:], radial_densities['on_target'], label="on target")
plt.plot(rs[1:], radial_densities['off_target'], label="off target")
plt.xlabel("Radial distance from Dot Center (Pixels)")
plt.ylabel("Average Dot Density Across all Hybs/Channels (Dots/Pixel)")
plt.legend()
plt.savefig(snakemake.output[1])
