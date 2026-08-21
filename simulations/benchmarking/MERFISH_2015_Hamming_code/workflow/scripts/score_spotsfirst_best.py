import pandas as pd
import numpy as np

truth = pd.read_csv(snakemake.input["truth"])
cb = pd.read_csv(snakemake.input["codebook"])

trth_cnts = truth.groupby("gene").size()
total_trth_cnts = np.sum(trth_cnts)


def add_truth_to_sumstats(sum_stats_path, decoded_path):
    sumstats = pd.read_csv(sum_stats_path)
    decoded = pd.read_csv(decoded_path)

    trth_cnts_df = pd.DataFrame({"truth": trth_cnts})
    if len(decoded) == 0:
        fp = 0
        trth_cnts_df["decoded"] = 0
        fdr = 0
    else:
        dcd_cnts = decoded.groupby("gene").size()
        dcd_cnts_df = pd.DataFrame({"decoded": dcd_cnts})
        dcd_cnts_df["truth"] = trth_cnts
        dcd_cnts_df.fillna(0, inplace=True)
        dcd_cnts_df["diff"] = dcd_cnts_df["truth"] - dcd_cnts_df["decoded"]

        fp = -np.sum(dcd_cnts_df["diff"].loc[dcd_cnts_df["diff"] <= 0])
        trth_cnts_df["decoded"] = dcd_cnts
        fdr = fp / np.sum(dcd_cnts)

    trth_cnts_df.fillna(0, inplace=True)
    trth_cnts_df["diff"] = trth_cnts_df["truth"] - trth_cnts_df["decoded"]

    not_decoded_true = np.sum(trth_cnts_df["diff"].loc[trth_cnts_df["diff"] >= 0])
    eff = (total_trth_cnts - not_decoded_true) / total_trth_cnts

    sumstats["eff_truth"] = eff
    sumstats["fdr_truth"] = fdr
    sumstats["fp"] = fp
    sumstats["not_decoded_true"] = not_decoded_true
    sumstats["total_true_cnts"] = total_trth_cnts
    sumstats["nbarcodes"] = snakemake.wildcards["nbarcodes"]
    sumstats["rep"] = snakemake.wildcards["rep"]
    sumstats["code"] = snakemake.wildcards["code"]
    sumstats["ncodewords"] = len(cb)
    sumstats["pdrop"] = snakemake.wildcards["pdrop"]

    return sumstats


scores = []
for sum_stats_path, decoded_path in zip(snakemake.input["sum_stats"], snakemake.input["decoded"]):
    scores.append(add_truth_to_sumstats(sum_stats_path, decoded_path))

all_scores = pd.concat(scores, ignore_index=True)

no_fd = all_scores[all_scores["fp"] == 0]

if len(no_fd) > 0:
    best = no_fd.loc[[no_fd["eff_truth"].idxmax()]]
    sensitivity = best["eff_truth"].iloc[0]
    fdr = best["fdr_truth"].iloc[0]
    n_true = best["total_true_cnts"].iloc[0]
    n_decoded = best["n_barcodes"].iloc[0]
else:
    sensitivity = 0.0
    fdr = np.nan
    n_true = all_scores["total_true_cnts"].iloc[0] if len(all_scores) > 0 else np.nan
    n_decoded = 0

result = pd.DataFrame({
    "method": [snakemake.params["method"]],
    "code": [snakemake.wildcards["code"]],
    "nbarcodes": [snakemake.wildcards["nbarcodes"]],
    "rstdv": [snakemake.wildcards["rstdv"]],
    "rep": [snakemake.wildcards["rep"]],
    "sensitivity": [sensitivity],
    "fdr": [fdr],
    "n_true": [n_true],
    "n_decoded": [n_decoded],
})

result.to_csv(snakemake.output[0], index=False)
