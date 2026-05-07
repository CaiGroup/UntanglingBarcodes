import pandas as pd
import matplotlib.pyplot as plt
import os
import skimage
import numpy as np
from scipy import stats
import scipy
import networkx as nx

def rm_dupes(gene_df, r):
    gene_df.reset_index(drop=True, inplace=True)
    gene_tree = scipy.spatial.KDTree(gene_df.loc[:, ['x', 'y']])
    pairs = gene_tree.query_pairs(r)
    g = nx.Graph()
    g.add_edges_from(pairs)
    ccs = nx.connected_components(g)
    duplicates = []
    for ds in ccs:
        ds = list(ds)
        largest_lambda_over_thresh = gene_df.iloc[ds, :].loc[:,'largest_lambda_over_thresh']
        ds.pop(np.argmax(largest_lambda_over_thresh))
        duplicates.extend(ds)
    drop_indices = gene_df.index[duplicates]
    return gene_df.drop(drop_indices) 

def plot_est_fdr_vs_sens(cp_nnc, cp_wnc, num_ge_cws, num_nc_cws, smFISH, smFISH_genes, rdup=1.6):
    cp_nnc_rmdp = pd.concat([rm_dupes(g, rdup) for _, g in cp_nnc.groupby('gene_number')])
    cp_wnc_rmdp = pd.concat([rm_dupes(g, rdup) for _, g in cp_wnc.groupby('gene_number')])

    lambdas = []
    sensitivities = []
    eFDRs = []
    ncs = []
    for lmbda in -np.sort(-np.unique(cp_nnc_rmdp.largest_lambda_over_thresh)):
        dcd_lmbda = cp_nnc_rmdp.loc[cp_nnc_rmdp.largest_lambda_over_thresh >= lmbda, :]
        dcd_gene_cnts = dcd_lmbda.groupby("gene").size()
        dcd_gene_cnts.index = [g.lower() for g in dcd_gene_cnts.index]
        ndcd = np.sum(np.array(dcd_gene_cnts))
        dcd_gene_cnts = dcd_gene_cnts.loc[[g in smFISH_genes for g in dcd_gene_cnts.index]]
        cp_wnc_rmdpdcd = cp_wnc_rmdp.loc[cp_wnc_rmdp.largest_lambda_over_thresh >= lmbda, :]
        lambda_nnc = np.sum(cp_wnc_rmdpdcd.gene == "negative_control")

        #print("len(dcd_gene_cnts): ", len(dcd_gene_cnts))
        if len(dcd_gene_cnts) > 0:
            results = pd.DataFrame({'seqFISH_ave_counts': dcd_gene_cnts, 'smFISH_ave_counts': smFISH}) #, index=dcd_gene_cnts.index)
            results.fillna(0, inplace=True)
            lm_synd_pub = stats.linregress(results['smFISH_ave_counts'], results['seqFISH_ave_counts'])
            sensitivities.append(lm_synd_pub.slope)
            eFDRs.append(num_ge_cws * (lambda_nnc / num_nc_cws) / ndcd)
        else:
            ("nothing!")
            sensitivities.append(0.0)
            eFDRs.append(0.0)

        lambdas.append(lmbda)
        ncs.append(lambda_nnc)

    sum_df = pd.DataFrame({'lambda': lambdas, 'sensitivity': sensitivities, 'ncs': ncs, 'est_fdr': eFDRs})
    return sum_df

def best_curve(pnts, res_row):
    return not any((res_row['est_fdr'] >= pnts['est_fdr']) & (res_row['sensitivity'] < pnts['sensitivity'])) 


results_ch488 = pd.read_csv(snakemake.input['bc488'])
cb488 = pd.read_csv(snakemake.input['cb488'])
n_ge_cws_488 = len(cb488)
n_nc_cws_488 = 8000 - n_ge_cws_488

results_ch561 = pd.read_csv(snakemake.input['bc561'])
cb561 = pd.read_csv(snakemake.input['cb561'])
n_ge_cws_561 = len(cb561)
n_nc_cws_561 = 8000 - n_ge_cws_561

results_ch643 = pd.read_csv(snakemake.input['bc643'])
cb643 = pd.read_csv(snakemake.input['cb643'])
n_ge_cws_643 = len(cb643)
n_nc_cws_643 = 8000 - n_ge_cws_643

results_ch488 = results_ch488.loc[results_ch488.sensitivity > 0, ['itga7', 'ksr1', 'vangl2',
       'tomm34', 'cd2ap', 'mlx', 'smarca4', 'nfkb2', 'crk', 'ccnd1', 'col4a1',
       'fam120a', 'trap1', 'slc25a1', 'ppp2ca', 'pdia3', 'eif4b', 'bgn',
       'lambda', 'sensitivity', 'ndcd', 'ncs', 'est_fdr', 'cell', 'beta_threshold']].reset_index()

results_ch561 = results_ch561.loc[results_ch561.sensitivity > 0, ['itga4', 'egfr', 'cbl',
       'ppp2r1b', 'ets2', 'axin1', 'sh3kbp1', 'itga5', 'ddx6', 'fxr2',
       'poldip3', 'sec61a1', 'itgb5', 'pdia4', 'eif3b', 'itgb1', 'serpinh1',
       'psat1', 'lambda', 'sensitivity', 'ndcd', 'ncs', 'est_fdr', 'cell',
       'beta_threshold']].reset_index()

results_ch643 = results_ch643.loc[results_ch643.sensitivity > 0, ['flt1', 'mxd1', 'cblb',
       'rabif', 'vangl1', 'tomm40l', 'dvl3', 'smarce1', 'dvl1', 'dvl2',
       'pdgfrb', 'smarcb1', 'trim8', 'cdkn1a', 'tomm40', 'myc', 'igfbp4',
       'ctnna1', 'vat1', 'zfp36l1', 'trim28', 'cyr61', 'fbln2', 'col1a1',
       'lambda', 'sensitivity', 'ndcd', 'ncs', 'est_fdr', 'cell', 'beta_threshold']].reset_index()

results_ch488["pos_cell"] = [c.split("_lf")[0] for c in results_ch488.cell]
results_ch561["pos_cell"] = [c.split("_lf")[0] for c in results_ch561.cell]
results_ch643["pos_cell"] = [c.split("_lf")[0] for c in results_ch643.cell]

def get_cell_closest_est_fdr(_cell_best_results, target_eFDR, ret_cols = "all"):
    cell_best_results = _cell_best_results.reset_index()
    isbest = [best_curve(cell_best_results, r) for i, r in cell_best_results.iterrows()]
    cell_best_results =  cell_best_results.loc[isbest]
    if ret_cols == "all":
        return cell_best_results.iloc[np.argmin(np.abs(cell_best_results.est_fdr - target_eFDR)), :]
    else:
        return cell_best_results.iloc[np.argmin(np.abs(cell_best_results.est_fdr - target_eFDR)), :][ret_cols]
    
smFISH_results = pd.read_csv(snakemake.input.smFISH_ref)
ncells_smFISH = len(smFISH_results)
smFISH_ave_counts = np.sum(smFISH_results, axis=0)/ncells_smFISH
smFISH_ave_counts.index.rename("gene", inplace=True)
smFISH_ave_counts.index = [g.lower() for g in smFISH_ave_counts.index]


def get_all_cell_fdr_sens(smFISH, results_ch488, results_ch561, results_ch643, target_fdr):
    def get_cefdr_488(c):
        return get_cell_closest_est_fdr(c, target_fdr,['itga7', 'ksr1', 'vangl2',
           'tomm34', 'cd2ap', 'mlx', 'smarca4', 'nfkb2', 'crk', 'ccnd1', 'col4a1',
           'fam120a', 'trap1', 'slc25a1', 'ppp2ca', 'pdia3', 'eif4b', 'bgn', 'ndcd', 'ncs'])
    
    def get_cefdr_561(c):
        return get_cell_closest_est_fdr(c, target_fdr,['itga4', 'egfr', 'cbl',
           'ppp2r1b', 'ets2', 'axin1', 'sh3kbp1', 'itga5', 'ddx6', 'fxr2',
           'poldip3', 'sec61a1', 'itgb5', 'pdia4', 'eif3b', 'itgb1', 'serpinh1',
           'psat1', 'ndcd', 'ncs'])
    
    def get_cefdr_643(c):
        return get_cell_closest_est_fdr(c, target_fdr,['flt1', 'mxd1', 'cblb',
           'rabif', 'vangl1', 'tomm40l', 'dvl3', 'smarce1', 'dvl1', 'dvl2',
           'pdgfrb', 'smarcb1', 'trim8', 'cdkn1a', 'tomm40', 'myc', 'igfbp4',
           'ctnna1', 'vat1', 'zfp36l1', 'trim28', 'cyr61', 'fbln2', 'col1a1', 'ndcd', 'ncs'])
        
    cr488 = results_ch488.groupby("pos_cell").apply(get_cefdr_488)
    crm488 = cr488.mean()
    
    cr561 = results_ch561.groupby("pos_cell").apply(get_cefdr_561)
    crm561 = cr561.mean()
    
    cr643 = results_ch643.groupby("pos_cell").apply(get_cefdr_643)
    crm643 = cr643.mean()

    est_fdr = (cr488.ncs.sum() + cr561.ncs.sum() + cr643.ncs.sum()) * (10000/14000)/(cr488.ndcd.sum() + cr561.ndcd.sum() + cr643.ndcd.sum())
    combined = pd.DataFrame({"smFISH" : smFISH_ave_counts, "seqFISH" : pd.concat([crm488, crm561, crm643])[smFISH_ave_counts.index]})

    lm_synd_pub = stats.linregress(combined.smFISH, combined.seqFISH)

    return est_fdr, lm_synd_pub.slope


def mapable_get_all_cell_fdr_sens(target_fdr):
    return get_all_cell_fdr_sens(smFISH_ave_counts, results_ch488, results_ch561, results_ch643, target_fdr)

fdr_sens = np.array(list(map(mapable_get_all_cell_fdr_sens, np.linspace(0.01, 0.2, 20))))

fdr_sens_curve = pd.DataFrame(fdr_sens, columns=["Estimated FDR", "Sensitivity relative to smFISH"])

plt.plot(fdr_sens_curve['Estimated FDR'], fdr_sens_curve['Sensitivity relative to smFISH'], label="Compressed Sensing Decoding")
plt.xlabel("Estimated FDR")
plt.ylabel("Sensitivity Relative to smFISH")
plt.savefig(snakemake.output[0])

fdr_sens_curve.to_csv(snakemake.output[1])

# Now plot sample regression showing the sensitivty calculation for estimated FDR = 0.05

target_fdr = 0.05

def get_cefdr_488(c):
    return get_cell_closest_est_fdr(c, target_fdr,['itga7', 'ksr1', 'vangl2',
       'tomm34', 'cd2ap', 'mlx', 'smarca4', 'nfkb2', 'crk', 'ccnd1', 'col4a1',
       'fam120a', 'trap1', 'slc25a1', 'ppp2ca', 'pdia3', 'eif4b', 'bgn', 'ndcd', 'ncs'])

def get_cefdr_561(c):
    return get_cell_closest_est_fdr(c, target_fdr,['itga4', 'egfr', 'cbl',
       'ppp2r1b', 'ets2', 'axin1', 'sh3kbp1', 'itga5', 'ddx6', 'fxr2',
       'poldip3', 'sec61a1', 'itgb5', 'pdia4', 'eif3b', 'itgb1', 'serpinh1',
       'psat1', 'ndcd', 'ncs'])

def get_cefdr_643(c):
    return get_cell_closest_est_fdr(c, target_fdr,['flt1', 'mxd1', 'cblb',
       'rabif', 'vangl1', 'tomm40l', 'dvl3', 'smarce1', 'dvl1', 'dvl2',
       'pdgfrb', 'smarcb1', 'trim8', 'cdkn1a', 'tomm40', 'myc', 'igfbp4',
       'ctnna1', 'vat1', 'zfp36l1', 'trim28', 'cyr61', 'fbln2', 'col1a1', 'ndcd', 'ncs'])
    
cr488 = results_ch488.groupby("pos_cell").apply(get_cefdr_488)
crm488 = cr488.mean()

cr561 = results_ch561.groupby("pos_cell").apply(get_cefdr_561)
crm561 = cr561.mean()

cr643 = results_ch643.groupby("pos_cell").apply(get_cefdr_643)
crm643 = cr643.mean()

estimated_FDR = (cr488.ncs.sum() + cr561.ncs.sum() + cr643.ncs.sum()) * (10000/14000)/(cr488.ndcd.sum() + cr561.ndcd.sum() + cr643.ndcd.sum())

combined = pd.DataFrame({"smFISH" : smFISH_ave_counts, "seqFISH" : pd.concat([crm488, crm561, crm643])[smFISH_ave_counts.index]})

lm_synd_pub = stats.linregress(combined.smFISH, combined.seqFISH)

plt.plot(combined.smFISH, combined.seqFISH, '.', label = "Mean transcripts per cell")
xs = np.array([np.min(combined.smFISH), np.max(combined.smFISH)])
plt.plot(xs, xs * lm_synd_pub.slope + lm_synd_pub.intercept, label = 'slope (sensivity) = ' + str(np.round(lm_synd_pub.slope, 2)))
plt.xlabel("smFISH")
plt.ylabel("seqFISH compressed sensing (estimated FDR = 0.05)")
plt.legend()
plt.savefig(snakemake.output[2])
