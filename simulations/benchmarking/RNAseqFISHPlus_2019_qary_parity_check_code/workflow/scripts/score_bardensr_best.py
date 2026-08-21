import pandas as pd
import numpy as np

truth = pd.read_csv(snakemake.input["truth"])
cb = pd.read_csv(snakemake.input["codebook"])
gene_col = cb.columns[0]

trth_cnts = truth.groupby(gene_col).size()
total_trth_cnts = np.sum(trth_cnts)


def score_decoded_file(decoded_path):
    decoded = pd.read_csv(decoded_path)

    if "gene" not in decoded.columns and "gene_number" in decoded.columns:
        decoded = decoded.rename(columns={"gene_number": "gene"})

    trth_cnts_df = pd.DataFrame({"truth": trth_cnts})

    if len(decoded) == 0:
        fd = 0
        trth_cnts_df["decoded"] = 0
        fdr = 0
        not_decoded_true = total_trth_cnts
    else:
        dcd_cnts = decoded.groupby("gene").size()
        dcd_cnts_df = pd.DataFrame({"decoded": dcd_cnts})
        dcd_cnts_df["truth"] = trth_cnts
        dcd_cnts_df.fillna(0, inplace=True)
        dcd_cnts_df["diff"] = dcd_cnts_df["truth"] - dcd_cnts_df["decoded"]

        fd = -np.sum(dcd_cnts_df["diff"].loc[dcd_cnts_df["diff"] <= 0])
        trth_cnts_df["decoded"] = dcd_cnts
        fdr = fd / np.sum(dcd_cnts) if np.sum(dcd_cnts) > 0 else 0

        trth_cnts_df.fillna(0, inplace=True)
        trth_cnts_df["diff"] = trth_cnts_df["truth"] - trth_cnts_df["decoded"]
        not_decoded_true = np.sum(trth_cnts_df["diff"].loc[trth_cnts_df["diff"] >= 0])

    sensitivity = (total_trth_cnts - not_decoded_true) / total_trth_cnts if total_trth_cnts > 0 else 0

    return pd.DataFrame({
        "method": [snakemake.params["method"]],
        "code": [snakemake.wildcards["code"]],
        "nbarcodes": [snakemake.wildcards["nbarcodes"]],
        "rstdv": [snakemake.wildcards["rstdv"]],
        "rep": [snakemake.wildcards["rep"]],
        "sensitivity": [sensitivity],
        "fdr": [fdr],
        "n_true": [int(total_trth_cnts)],
        "n_decoded": [len(decoded)],
        "n_fd": [int(fd)],
        "not_decoded_true": [int(not_decoded_true)]
    })


scores = []
for f, (l1, peak) in zip(snakemake.input["decoded"], snakemake.params["combos"]):
    df = score_decoded_file(f)
    df["param_l1"] = l1
    df["param_peak"] = peak
    scores.append(df)

all_scores = pd.concat(scores, ignore_index=True)

no_fd = all_scores[all_scores["n_fd"] == 0]

if len(no_fd) > 0:
    best = no_fd.loc[[no_fd["sensitivity"].idxmax()]].copy()
else:
    best = pd.DataFrame({
        "method": [snakemake.params["method"]],
        "code": [snakemake.wildcards["code"]],
        "nbarcodes": [snakemake.wildcards["nbarcodes"]],
        "rstdv": [snakemake.wildcards["rstdv"]],
        "rep": [snakemake.wildcards["rep"]],
        "sensitivity": [0.0],
        "fdr": [np.nan],
        "n_true": [all_scores["n_true"].iloc[0] if len(all_scores) > 0 else np.nan],
        "n_decoded": [0],
        "n_fd": [np.nan],
        "not_decoded_true": [np.nan],
        "param_l1": [np.nan],
        "param_peak": [np.nan],
    })

best.to_csv(snakemake.output[0], index=False)
