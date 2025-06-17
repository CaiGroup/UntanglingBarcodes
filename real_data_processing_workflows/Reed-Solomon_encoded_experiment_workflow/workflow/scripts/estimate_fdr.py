import pandas as pd
import numpy as np

ss_w_neg_ctrl = pd.read_csv(snakemake.input[0])
ss_wo_neg_ctrl = pd.read_csv(snakemake.input[1])
cb = pd.read_csv(snakemake.input[2])
if len(ss_w_neg_ctrl) == 0:
    print("len zeros")
    ss_wo_neg_ctrl.to_csv(snakemake.output[0])
else:
    print("ss_w_neg_ctrl")
    print(ss_w_neg_ctrl)

    print("ss_wo_neg_ctrl")
    print(ss_wo_neg_ctrl)

    n_neg_ctrl_cws = np.sum(cb.gene == "negative_control")
    n_gene_encoding_cws = len(cb) - n_neg_ctrl_cws

    n_neg_ctrl_barcodes = ss_w_neg_ctrl.negative_control_barcodes[0]

    est_fdr_per_cw = n_neg_ctrl_barcodes/n_neg_ctrl_cws
    est_fds = est_fdr_per_cw * n_gene_encoding_cws
    est_fdr = est_fds/ss_wo_neg_ctrl.gene_encoding_barcodes[0]

    ss_wo_neg_ctrl.loc[0,"est_fp_rate"] = est_fdr
    ss_wo_neg_ctrl.loc[0,"negative_control_barcodes"] = n_neg_ctrl_barcodes

    print("saving...")
    print(ss_wo_neg_ctrl)
    ss_wo_neg_ctrl.to_csv(snakemake.output[0])
