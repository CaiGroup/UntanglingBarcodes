import pandas as pd
import numpy as np
import re
import os

scores = []
for f in snakemake.input:
    df = pd.read_csv(f)
    m = re.search(r'_k1([\d.]+)\.csv$', os.path.basename(f))
    if m:
        df["param_k1"] = float(m.group(1))
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
        "param_k1": [np.nan],
    })

best.to_csv(snakemake.output[0], index=False)
