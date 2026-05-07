import sys
import pandas as pd
import matplotlib.pyplot as plt
import os
import skimage
import numpy as np
from scipy import stats
import scipy
import networkx as nx

def best_curve(pnts, res_row):
    return not any((res_row['est_fdr'] >= pnts['est_fdr']) & (res_row['sensitivity'] < pnts['sensitivity'])) 

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

def plot_est_fdr_vs_sens(cp_nnc, cp_wnc, num_ge_cws, num_nc_cws, smFISH, rdup=1.6):
    print("ncpaths nnc before rm dps:", len(cp_nnc))
    cp_nnc_rmdp = pd.concat([rm_dupes(g, rdup) for _, g in cp_nnc.groupby('gene_number')])
    cp_wnc_rmdp = pd.concat([rm_dupes(g, rdup) for _, g in cp_wnc.groupby('gene_number')])

    lambdas = []
    sensitivities = []
    eFDRs = []
    ncs = []
    gene_cnts = []
    ndcds = []
    for lmbda in -np.sort(-np.unique(cp_nnc_rmdp.largest_lambda_over_thresh)):
        dcd_lmbda = cp_nnc_rmdp.loc[cp_nnc_rmdp.largest_lambda_over_thresh >= lmbda, :]
        dcd_gene_cnts = dcd_lmbda.groupby("gene").size()
        dcd_gene_cnts.index = [g.lower() for g in dcd_gene_cnts.index]
        ndcd = np.sum(np.array(dcd_gene_cnts))
        dcd_gene_cnts = dcd_gene_cnts.loc[[g in smFISH.index for g in dcd_gene_cnts.index]]
        cp_wnc_rmdpdcd = cp_wnc_rmdp.loc[cp_wnc_rmdp.largest_lambda_over_thresh >= lmbda, :]
        lambda_nnc = np.sum(cp_wnc_rmdpdcd.gene == "negative_control")

        #print("len(dcd_gene_cnts): ", len(dcd_gene_cnts))
        if len(dcd_gene_cnts) > 0:
            results = pd.DataFrame({'seqFISH_ave_counts': dcd_gene_cnts, 'smFISH_ave_counts': smFISH}) #, index=dcd_gene_cnts.index)
            results.fillna(0, inplace=True)
            lm_synd_pub = stats.linregress(results['smFISH_ave_counts'], results['seqFISH_ave_counts'])
            sensitivities.append(lm_synd_pub.slope)
            eFDRs.append(num_ge_cws * (lambda_nnc / num_nc_cws) / ndcd)
            gene_cnts.append(pd.DataFrame(results['seqFISH_ave_counts']).T)
            #print("results['seqFISH_ave_counts'].T")
            #print(pd.DataFrame(results['seqFISH_ave_counts']).T)
        else:
            ("nothing!")
            #print("pd.DataFrame(pd.Series(0, smFISH.index)).T")
            #print(pd.DataFrame(pd.Series(0, smFISH.index)).T)
            gene_cnts.append(pd.DataFrame(pd.Series(0, smFISH.index)).T)
            sensitivities.append(0.0)
            eFDRs.append(0.0)

        lambdas.append(lmbda)
        ncs.append(lambda_nnc)
        ndcds.append(ndcd)

    sum_df_gene_cnts = pd.concat(gene_cnts)
    sum_df_gene_cnts.reset_index(drop=True, inplace=True)
    #print('sum_df_gene_cnts')
    #print(sum_df_gene_cnts)
    #print( pd.DataFrame({'lambda': lambdas, 'sensitivity': sensitivities, 'ncs': ncs, 'est_fdr': eFDRs}))
    sum_df = pd.concat([sum_df_gene_cnts, pd.DataFrame({'lambda': lambdas, 'sensitivity': sensitivities, 'ndcd' : ndcds, 'ncs': ncs, 'est_fdr': eFDRs})], axis=1)
    #print('sum_df')
    #print(sum_df)
    return sum_df

def main():
    nnc_files = snakemake.input.nnc_files
    wnc_files = snakemake.input.wnc_files
    lfs = snakemake.params.lfs
    search_radii = snakemake.params.search_radii
    smFISH_results = pd.read_csv(snakemake.input.smFISH_ref)
    #smFISH_results.set_index("cell id", inplace=True)
    
    ncells_smFISH = len(smFISH_results)
    smFISH_ave_counts = np.sum(smFISH_results, axis=0)/ncells_smFISH
    smFISH_ave_counts.index.rename("gene", inplace=True)

    smFISH_ave_counts.index = [g.lower() for g in smFISH_ave_counts.index]


    codebook = pd.read_csv(snakemake.input.codebook)

    codebook_genes = [g.lower() for g in codebook['gene']]

    smFISH_ave_counts = smFISH_ave_counts.loc[[g in codebook_genes for g in smFISH_ave_counts.index]]

    num_ge_cws = np.sum(codebook["gene"] != "negative_control")
    num_nc_cws = np.sum(codebook["gene"] == "negative_control")

    nnc_betas = snakemake.input.nnc_beta_paths 
    wnc_betas = snakemake.input.wnc_beta_paths 
    beta_thresholds = snakemake.params.beta_thresholds

    path_summary_files = snakemake.input.paths_summaries


    # Organize files by search radius
    files_by_sr = {str(sr): [] for sr in search_radii}
    for nnc, wnc, b_nnc, b_wnc, lpsum in zip(nnc_files, wnc_files, nnc_betas, wnc_betas, path_summary_files):
    # Extract sr value from filename
        sr_val = None
        for sr in search_radii:
            if f"_sr{sr}" in nnc:
                sr_val = str(sr)
                break
        if sr_val is not None:
            files_by_sr[sr_val].append((nnc, wnc, b_nnc, b_wnc, lpsum))

    n_sr = len(search_radii)
    n_bta_thr = len(beta_thresholds)
    fig, axes = plt.subplots(n_bta_thr, n_sr, figsize=(7 * n_sr, 6*n_bta_thr), squeeze=False)
    sum_dfs = []

    for idx, sr in enumerate(search_radii):
        sr_str = str(sr)
        for lf_val in lfs:
            # Find files for this lf and sr
            nnc_path = next((f[0] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[0]), None)
            wnc_path = next((f[1] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[1]), None)
            nnc_b_path = next((f[2] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[2]), None)
            wnc_b_path = next((f[3] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[3]), None)
            lp_summary_path = next((f[4] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[4]), None)
            if nnc_path is None or wnc_path is None or nnc_b_path is None or wnc_b_path is None:
                continue
            if not (os.path.exists(nnc_path) and os.path.exists(wnc_path) and os.path.exists(nnc_b_path) and os.path.exists(wnc_b_path)):
                continue
            cp_nnc = pd.read_csv(nnc_path)
            cp_wnc = pd.read_csv(wnc_path)
            b_nnc = skimage.io.imread(nnc_b_path)
            b_wnc = skimage.io.imread(wnc_b_path)
            lp_summary = pd.read_csv(lp_summary_path)

            cell_str = nnc_path.split("nnc_")[1].split(".csv")[0]
            
            for jdx, bta_thr in enumerate(beta_thresholds):
                ax = axes[jdx, idx]

                b_thresh_nnc = np.any(b_nnc >= bta_thr, axis=1)
                b_thresh_wnc = np.any(b_wnc >= bta_thr, axis=1)

                cp_nnc_bt = cp_nnc.loc[b_thresh_nnc,:]
                cp_wnc_bt = cp_wnc.loc[b_thresh_wnc,:]

                cp_nnc_bt.loc[:, 'largest_lambda_over_thresh'] = np.array(lp_summary['lambda'][(b_nnc[b_thresh_nnc,:] >= bta_thr).argmax(axis=1)])
                cp_wnc_bt.loc[:, 'largest_lambda_over_thresh'] = np.array(lp_summary['lambda'][(b_wnc[b_thresh_wnc,:] >= bta_thr).argmax(axis=1)])

                if len(cp_nnc_bt) > 0 and len(cp_wnc_bt) > 0:
                
                    sum_df = plot_est_fdr_vs_sens(cp_nnc_bt, cp_wnc_bt, num_ge_cws, num_nc_cws, smFISH_ave_counts, rdup=float(sr))
                    if lf_val > 2.5:
                        ax.plot(sum_df.loc[5:, 'est_fdr'], sum_df.loc[5:, 'sensitivity'],'--', marker='.', label=f'lf={lf_val}')
                    else:
                        ax.plot(sum_df.loc[5:, 'est_fdr'], sum_df.loc[5:, 'sensitivity'], marker='.', label=f'lf={lf_val}')
                    ax.set_xlabel('Estimated FDR')
                    ax.set_ylabel('Sensitivity relative to smFISH')
                    ax.set_title(f'Search Radius = {sr}, Beta_thresh = {bta_thr}')
                    ax.legend()
                    sum_df['cell'] = cell_str
                    sum_df['beta_threshold'] = bta_thr
                    sum_dfs.append(sum_df)


    print("saving..")
    plt.tight_layout()
    plt.savefig(snakemake.output[0])
    plt.close()
    print("saved")
    combined_sum_df = pd.concat(sum_dfs)
    isbest = [best_curve(combined_sum_df, r) for i, r in combined_sum_df.iterrows()]

    best_sum_df = combined_sum_df.iloc[isbest, :]
    plt.plot(best_sum_df['est_fdr'], best_sum_df['sensitivity'])
    plt.xlabel('Estimated FDR')
    plt.ylabel('Sensitivity relative to smFISH')
    plt.savefig(snakemake.output[1])
    best_sum_df.to_csv(snakemake.output[2], index=False)
    
if __name__ == "__main__":
    main()