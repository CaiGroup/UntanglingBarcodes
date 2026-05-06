import pandas as pd
import numpy as np

sumstats = pd.read_csv(snakemake.input[0])
decoded = pd.read_csv(snakemake.input[1])
truth = pd.read_csv(snakemake.input[2])
cb = pd.read_csv(snakemake.input[3])

ngenes = np.sum(cb.gene != "negative_control")

dcd_cnts = decoded.groupby("gene_number").size()
trth_cnts = truth.groupby("gene_number").size()

total_trth_cnts = np.sum(trth_cnts)

dcd_cnts_df = pd.DataFrame({"decoded": dcd_cnts})
trth_cnts_df = pd.DataFrame({"truth": trth_cnts})
dcd_cnts_df= dcd_cnts_df.loc[dcd_cnts.index <= ngenes]

dcd_cnts_df["truth"] = trth_cnts
trth_cnts_df["decoded"] = dcd_cnts

dcd_cnts_df.fillna(0, inplace=True)
trth_cnts_df.fillna(0, inplace=True)

trth_cnts_df["diff"] = trth_cnts_df["truth"] - trth_cnts_df["decoded"]

not_decoded_true = np.sum(trth_cnts_df["diff"].loc[trth_cnts_df["diff"] >=0])
eff = (total_trth_cnts - not_decoded_true)/total_trth_cnts


dcd_cnts_df["diff"] = dcd_cnts_df["truth"] - dcd_cnts_df["decoded"]

fp = -np.sum(dcd_cnts_df["diff"].loc[dcd_cnts_df["diff"] <= 0])

fdr = fp/np.sum(dcd_cnts)

sumstats["eff_truth"] = eff
sumstats["fdr_truth"] = fdr
sumstats["fp"] = fp
sumstats["not_decoded_true"] = not_decoded_true
sumstats["total_true_cnts"] = total_trth_cnts
sumstats["nbarcodes"] = snakemake.wildcards['nbarcodes']
sumstats["nnonspec"] = snakemake.wildcards['nnonspec']
sumstats["ims"] = snakemake.wildcards["ims"]
sumstats["rep"] = snakemake.wildcards["rep"]

sumstats.to_csv(snakemake.output[0])