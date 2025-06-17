import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input[0])

filt_size = snakemake.params["filt_size"]
filt_prop = snakemake.params["filt_prop"]

def neg_cntrl_saturated_ccs(df):
    n_nctrls = np.sum(df.gene == "negative_control")
    remove_if_true = len(df) > 5 and n_nctrls/len(df) > 0.3
    #print("remove_if_true: ",  remove_if_true)
    return not remove_if_true

df_filtered = df.groupby("cc").filter(neg_cntrl_saturated_ccs)

df_filtered.to_csv(snakemake.output[0])

# calculate summary stats

n_barcodes = len(df_filtered)
n_neg_ctrl_barcodes = np.sum(df_filtered.gene == 'negative_control')

cb = pd.read_csv(snakemake.input[1])
n_codewords = len(cb)
n_ontargets = n_codewords - np.sum(cb.gene == "negative_control")

sum_stats = pd.DataFrame()
sum_stats.loc[:,"pos"] = [snakemake.wildcards['pos']]
sum_stats.loc[:, "ch"] = [snakemake.wildcards['ch']]
sum_stats.loc[:,"lvf"] = [snakemake.wildcards['lf']]
sum_stats.loc[:,"lwvf"] = [snakemake.wildcards['wf']]
sum_stats.loc[:,"svf"] = [snakemake.wildcards['sf']]
sum_stats.loc[:,"dr"] = [snakemake.wildcards['dr']]
sum_stats.loc[:,"n_barcodes"] = [n_barcodes]
sum_stats.loc[:,"negative_control_barcodes"] = [n_neg_ctrl_barcodes]
sum_stats.loc[:,"gene_encoding_barcodes"] = [n_barcodes - n_neg_ctrl_barcodes]
sum_stats.loc[:,"est_fp_rate"] = [(n_neg_ctrl_barcodes*n_ontargets/(n_codewords-n_ontargets))/(n_barcodes-n_neg_ctrl_barcodes)]

sum_stats.to_csv(snakemake.output[1])
