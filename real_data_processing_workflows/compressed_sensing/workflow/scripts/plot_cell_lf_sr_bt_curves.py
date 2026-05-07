import sys
import pandas as pd
import matplotlib.pyplot as plt
import os
import skimage
import numpy as np

def rm_dupes(gene_df, r):
    import scipy
    import networkx as nx
    gene_tree = scipy.spatial.KDTree(gene_df.loc[:, ['x', 'y']])
    pairs = gene_tree.query_pairs(r)
    g = nx.Graph()
    g.add_edges_from(pairs)
    ccs = nx.connected_components(g)
    duplicates = []
    for ds in ccs:
        ds = list(ds)
        ds.sort()
        for i in ds[1:]:
            duplicates.append(i)
    drop_indices = gene_df.index[duplicates]
    return gene_df.drop(drop_indices) 

def plot_est_fdr_vs_ge(cp_nnc, cp_wnc, rdup=1.6):
    cp_nnc_rmdp = pd.concat([rm_dupes(g, rdup) for _, g in cp_nnc.groupby('gene_number')])
    cp_wnc_rmdp = pd.concat([rm_dupes(g, rdup) for _, g in cp_wnc.groupby('gene_number')])
    marge_bcs_nnc = cp_nnc_rmdp.groupby('largest_lambda_over_thresh').size()
    cp_wc_rmdp = cp_wnc_rmdp.loc[cp_wnc_rmdp.gene == "negative_control", :]
    marge_ncs = cp_wc_rmdp.groupby('largest_lambda_over_thresh').size()
    marge_bcs_nnc.name = 'marge_bcs_nnc'
    sum_df = pd.DataFrame(marge_bcs_nnc)
    sum_df['marge_ncs'] = marge_ncs
    sum_df.fillna(0, inplace=True)
    sum_df.sort_index(ascending=False, inplace=True)
    sum_df["ge"] = sum_df.marge_bcs_nnc.cumsum()
    sum_df["nc"] = sum_df.marge_ncs.cumsum()
    sum_df['est_fdr'] = 3333 * (sum_df['nc'] / 4667) / sum_df['ge']
    return sum_df

def main():
    nnc_files = snakemake.input.nnc_files
    wnc_files = snakemake.input.wnc_files
    lfs = snakemake.params.lfs
    search_radii = snakemake.params.search_radii
    output_path = snakemake.output[0]


    #nnc_files = [f"results/lassoed_cpaths/lassoed_cpaths_nnc_cell_1_lf{lf}_sr{sr}.csv" for lf in lfs for sr in search_radii]
    #wnc_files = [f"results/lassoed_cpaths/lassoed_cpaths_wnc_cell_1_lf{lf}_sr{sr}.csv" for lf in lfs for sr in search_radii]
    nnc_betas = snakemake.input.nnc_beta_paths #[f"results/lassoed_cpaths/betas_nnc_cell_1_lf{lf}_sr{sr}.png" for lf in lfs for sr in search_radii]
    wnc_betas = snakemake.input.wnc_beta_paths #[f"results/lassoed_cpaths/betas_wnc_cell_1_lf{lf}_sr{sr}.png" for lf in lfs for sr in search_radii]
    #lfs = snakemake.params.lfs
    #search_radii = snakemake.params.search_radii
    beta_thresholds = snakemake.params.beta_thresholds
    #output_path = snakemake.output[0]

    # Organize files by search radius
    files_by_sr = {str(sr): [] for sr in search_radii}
    for nnc, wnc, b_nnc, b_wnc in zip(nnc_files, wnc_files, nnc_betas, wnc_betas):
    # Extract sr value from filename
        sr_val = None
        for sr in search_radii:
            if f"_sr{sr}" in nnc:
                sr_val = str(sr)
                break
        if sr_val is not None:
            files_by_sr[sr_val].append((nnc, wnc, b_nnc, b_wnc))

    n_sr = len(search_radii)
    n_bta_thr = len(beta_thresholds)
    fig, axes = plt.subplots(n_bta_thr, n_sr, figsize=(7 * n_sr, 6*n_bta_thr), squeeze=False)

    for idx, sr in enumerate(search_radii):
        sr_str = str(sr)
        for lf_val in lfs:
            # Find files for this lf and sr
            nnc_path = next((f[0] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[0]), None)
            wnc_path = next((f[1] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[1]), None)
            nnc_b_path = next((f[2] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[2]), None)
            wnc_b_path = next((f[3] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[3]), None)
            if nnc_path is None or wnc_path is None or nnc_b_path is None or wnc_b_path is None:
                continue
            if not (os.path.exists(nnc_path) and os.path.exists(wnc_path) and os.path.exists(nnc_b_path) and os.path.exists(wnc_b_path)):
                continue
            cp_nnc = pd.read_csv(nnc_path)
            cp_wnc = pd.read_csv(wnc_path)
            b_nnc = skimage.io.imread(nnc_b_path)
            b_wnc = skimage.io.imread(wnc_b_path)
            
            
            for jdx, bta_thr in enumerate(beta_thresholds):
                ax = axes[jdx, idx]

                b_thresh_nnc = np.any(b_nnc > bta_thr, axis=1)
                b_thresh_wnc = np.any(b_wnc > bta_thr, axis=1)

                cp_nnc_bt = cp_nnc.loc[b_thresh_nnc,:]
                cp_wnc_bt = cp_wnc.loc[b_thresh_wnc,:]
                
                sum_df = plot_est_fdr_vs_ge(cp_nnc_bt, cp_wnc_bt, rdup=float(sr))
                if lf_val > 2.5:
                    ax.plot(sum_df.loc[5:, 'est_fdr'], sum_df.loc[5:, 'ge'],'--', marker='.', label=f'lf={lf_val}')
                else:
                    ax.plot(sum_df.loc[5:, 'est_fdr'], sum_df.loc[5:, 'ge'], marker='.', label=f'lf={lf_val}')
                ax.set_xlabel('Estimated FDR')
                ax.set_ylabel('Decoded Gene Encoding Barcodes')
                ax.set_title(f'Search Radius = {sr}, Beta_thresh = {bta_thr}')
                ax.legend()


    print("saving..")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print("saved")
if __name__ == "__main__":
    main()