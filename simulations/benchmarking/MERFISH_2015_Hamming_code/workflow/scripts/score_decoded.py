import pandas as pd
import numpy as np

decoded = pd.read_csv(snakemake.input["decoded"])
truth = pd.read_csv(snakemake.input["truth"])
cb = pd.read_csv(snakemake.input["codebook"])

gene_col = cb.columns[0]
n_codewords = len(cb)

if "gene" not in decoded.columns and "gene_number" in decoded.columns:
    decoded = decoded.rename(columns={"gene_number": "gene"})

trth_cnts = truth.groupby(gene_col).size()
total_trth_cnts = np.sum(trth_cnts)
trth_cnts_df = pd.DataFrame({"truth": trth_cnts})

if len(decoded) == 0:
    fd = 0
    trth_cnts_df["decoded"] = 0
    fdr = 0
    not_decoded_true = total_trth_cnts
else:
    decoded_gene_col = "gene"
    dcd_cnts = decoded.groupby(decoded_gene_col).size()
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

result = pd.DataFrame({
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

result.to_csv(snakemake.output[0], index=False)
