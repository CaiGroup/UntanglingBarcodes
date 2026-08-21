import pandas as pd
import numpy as np

# Score the CS L0 decoding by the ACTUAL support at each lambda on the L0 path, and report the best
# zero-false-discovery sensitivity found at any lambda. This differs from score_cs_ranked.py, which
# ranks each codepath by its first-entry lambda and sweeps a monotone threshold -- an assumption the
# L0 support violates (codepaths can leave and re-enter), which can hide a clean operating point.

perlam = pd.read_csv(snakemake.input["decoded_perlambda"])
truth = pd.read_csv(snakemake.input["truth"])
cb = pd.read_csv(snakemake.input["codebook"])

gene_col = cb.columns[0]
trth_cnts = truth.groupby(gene_col).size()
total_trth_cnts = int(np.sum(trth_cnts))

best_sensitivity = 0.0
best_lambda = np.nan
n_decoded_at_best = 0

if len(perlam) > 0 and not perlam["lambda"].isna().all():
    for lam, grp in perlam.groupby("lambda"):
        # decoded count per gene = number of codepaths selected for that gene at this lambda
        dcd_cnts = grp.groupby(gene_col).size()
        cnts = pd.DataFrame({"decoded": dcd_cnts, "truth": trth_cnts})
        cnts.fillna(0, inplace=True)
        cnts["diff"] = cnts["truth"] - cnts["decoded"]
        # false discoveries = over-decoded genes (incl. spurious genes not in truth)
        fd = int(-np.sum(cnts["diff"].loc[cnts["diff"] <= 0]))
        if fd == 0:
            not_decoded_true = int(np.sum(cnts["diff"].loc[cnts["diff"] >= 0]))
            sensitivity = (total_trth_cnts - not_decoded_true) / total_trth_cnts if total_trth_cnts > 0 else 0.0
            if sensitivity > best_sensitivity:
                best_sensitivity = sensitivity
                best_lambda = lam
                n_decoded_at_best = len(grp)

result = pd.DataFrame({
    "method": [snakemake.params["method"]],
    "code": [snakemake.wildcards["code"]],
    "nbarcodes": [snakemake.wildcards["nbarcodes"]],
    "rstdv": [snakemake.wildcards["rstdv"]],
    "rep": [snakemake.wildcards["rep"]],
    "sensitivity": [best_sensitivity],
    "best_lambda": [best_lambda],
    "n_true": [total_trth_cnts],
    "n_decoded": [n_decoded_at_best],
})

result.to_csv(snakemake.output[0], index=False)
