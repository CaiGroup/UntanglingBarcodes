import pandas as pd
import numpy as np

decoded_ranked = pd.read_csv(snakemake.input["decoded_ranked"])
truth = pd.read_csv(snakemake.input["truth"])
cb = pd.read_csv(snakemake.input["codebook"])

gene_col = cb.columns[0]

trth_cnts = truth.groupby(gene_col).size()
total_trth_cnts = int(np.sum(trth_cnts))

if len(decoded_ranked) == 0 or decoded_ranked["lambda"].isna().all():
    sensitivity = 0.0
    best_lambda = np.nan
    n_decoded_at_best = 0
else:
    lambdas = sorted(decoded_ranked["lambda"].dropna().unique())

    best_lambda = np.nan
    sensitivity = 0.0
    n_decoded_at_best = 0

    for lam in lambdas:
        subset = decoded_ranked[decoded_ranked["lambda"] >= lam]

        dcd_cnts = subset.groupby(gene_col).size()
        dcd_cnts_df = pd.DataFrame({"decoded": dcd_cnts, "truth": trth_cnts})
        dcd_cnts_df.fillna(0, inplace=True)
        dcd_cnts_df["diff"] = dcd_cnts_df["truth"] - dcd_cnts_df["decoded"]

        fd = int(-np.sum(dcd_cnts_df["diff"].loc[dcd_cnts_df["diff"] <= 0]))

        if fd == 0:
            best_lambda = lam
            trth_cnts_df = pd.DataFrame({"truth": trth_cnts, "decoded": dcd_cnts})
            trth_cnts_df.fillna(0, inplace=True)
            trth_cnts_df["diff"] = trth_cnts_df["truth"] - trth_cnts_df["decoded"]
            not_decoded_true = int(np.sum(trth_cnts_df["diff"].loc[trth_cnts_df["diff"] >= 0]))
            sensitivity = (total_trth_cnts - not_decoded_true) / total_trth_cnts if total_trth_cnts > 0 else 0.0
            n_decoded_at_best = len(subset)
            break

result = pd.DataFrame({
    "method": [snakemake.params["method"]],
    "code": [snakemake.wildcards["code"]],
    "nbarcodes": [snakemake.wildcards["nbarcodes"]],
    "rstdv": [snakemake.wildcards["rstdv"]],
    "rep": [snakemake.wildcards["rep"]],
    "sensitivity": [sensitivity],
    "best_lambda": [best_lambda],
    "n_true": [total_trth_cnts],
    "n_decoded": [n_decoded_at_best],
})

result.to_csv(snakemake.output[0], index=False)
