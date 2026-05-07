import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import scipy
import networkx as nx
import skimage


def get_results_series_concat_best(lfs, search_radii, beta_thresholds):
    # --- Helper functions from notebook ---
    def rm_dupes(gene_df, r):
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

    def best_curve(pnts, res_row):
        return not any((res_row['est_fdr'] >= pnts['est_fdr']) & (res_row['ge'] < pnts['ge'])) 

    # --- File loading logic (requires snakemake and files to exist) ---
    nnc_files = snakemake.input["nnc_files"]
    nnc_betas = snakemake.input["nnc_betas"]
    wnc_files = snakemake.input["wnc_files"]
    wnc_betas = snakemake.input["wnc_betas"]
    
    print(nnc_files)
    print(nnc_betas)

    print(wnc_files)
    print(wnc_betas)


    # Organize files by search radius
    files_by_sr = {str(sr): [] for sr in search_radii}
    for nnc, wnc, b_nnc, b_wnc in zip(nnc_files, wnc_files, nnc_betas, wnc_betas): 
        sr_val = None
        for sr in search_radii:
            if f"_sr{sr}" in nnc:
                sr_val = str(sr)
                break
        if True:  # sr_val is not None:
            files_by_sr[sr_val].append((nnc, wnc, b_nnc, b_wnc))

    results_series = []
    print(files_by_sr)
    for idx, sr in enumerate(search_radii):
        sr_str = str(sr)
        for lf_val in lfs:
            print(f"_lf{lf_val}_sr{sr}")
            nnc_path = next((f[0] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[0]), None)
            wnc_path = next((f[1] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[1]), None)
            nnc_b_path = next((f[2] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[2]), None)
            wnc_b_path = next((f[3] for f in files_by_sr[sr_str] if f"_lf{lf_val}_sr{sr}" in f[3]), None)
            #nnc_path = next((f for f in files_by_sr[sr_str][0][0] if f"_lf{lf_val}_sr{sr}" in f), None)
            #wnc_path = next((f for f in files_by_sr[sr_str][0][1] if f"_lf{lf_val}_sr{sr}" in f), None)
            #nnc_b_path = next((f for f in files_by_sr[sr_str][0][2] if f"_lf{lf_val}_sr{sr}" in f), None)
            #wnc_b_path = next((f for f in files_by_sr[sr_str][0][3] if f"_lf{lf_val}_sr{sr}" in f), None)
            if nnc_path is None or wnc_path is None or nnc_b_path is None or wnc_b_path is None:
                print(nnc_path, wnc_path, nnc_b_path, wnc_b_path)
                print("all none")
                continue
            if not (os.path.exists(nnc_path) and os.path.exists(wnc_path) and os.path.exists(nnc_b_path) and os.path.exists(wnc_b_path)):
                print(nnc_path, wnc_path, nnc_b_path, wnc_b_path)
                print("paths don't exist")
                continue
            cp_nnc = pd.read_csv(nnc_path)
            cp_wnc = pd.read_csv(wnc_path)
            b_nnc = skimage.io.imread(nnc_b_path)
            b_wnc = skimage.io.imread(wnc_b_path)
            for jdx, bta_thr in enumerate(beta_thresholds):
                b_thresh_nnc = np.any(b_nnc > bta_thr, axis=1)
                b_thresh_wnc = np.any(b_wnc > bta_thr, axis=1)
                cp_nnc_bt = cp_nnc.loc[b_thresh_nnc, :]
                cp_wnc_bt = cp_wnc.loc[b_thresh_wnc, :]
                sum_df = plot_est_fdr_vs_ge(cp_nnc_bt, cp_wnc_bt, rdup=float(sr))
                sum_df['sr'] = sr
                sum_df['lf'] = lf_val
                sum_df['bta_thr'] = bta_thr
                results_series.append(sum_df)

    results_series_concat = pd.concat(results_series)
    isbest = [best_curve(results_series_concat, r) for i, r in results_series_concat.iterrows()]
    results_series_concat['best'] = isbest
    results_series_concat_best = results_series_concat.loc[results_series_concat.best, :]
    return results_series_concat_best


lfs = snakemake.params['lfs'] 
search_radii = snakemake.params['search_radii']
beta_thresholds = snakemake.params['beta_thresholds']

results_series_concat_best = get_results_series_concat_best(lfs, search_radii, beta_thresholds)

results_series_concat_best.to_csv(snakemake.output[0])
#dotsfirst_stats = pd.concat([pd.read_csv(fn) for fn in snakemake.input['dotsfirst_stats']])


#plt.plot(dotsfirst_stats['est_fp_rate'], dotsfirst_stats['gene_encoding_barcodes'], '.', label="dots first")
#plt.plot(results_series_concat_best.est_fdr, results_series_concat_best['ge'], '.')
#plt.xlabel('Estimated FDR')
#plt.ylabel('Number of decoded gene encoding barcodes')

#plt.tight_layout()
#plt.savefig(snakemake.output[0])
#plt.close()