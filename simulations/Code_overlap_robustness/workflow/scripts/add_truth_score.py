import pandas as pd
import numpy as np

sumstats = pd.read_csv(snakemake.input[0])
decoded = pd.read_csv(snakemake.input[1])
truth = pd.read_csv(snakemake.input[2])
cb = pd.read_csv(snakemake.input[3])

ngenes = len(cb) #np.sum(cb.gene != "negative_control")


trth_cnts = truth.groupby("gene").size()
total_trth_cnts = np.sum(trth_cnts)

trth_cnts_df = pd.DataFrame({"truth": trth_cnts})
if len(decoded) == 0:
    fp = 0
    trth_cnts_df["decoded"] = 0
    fdr = 0
else:
    dcd_cnts = decoded.groupby("gene").size()
    dcd_cnts_df = pd.DataFrame({"decoded": dcd_cnts})
    #print("ngenes: ", ngenes)
    #dcd_cnts_df= dcd_cnts_df.loc[dcd_cnts.index <= ngenes]
    dcd_cnts_df["truth"] = trth_cnts
    dcd_cnts_df.fillna(0, inplace=True)
    dcd_cnts_df["diff"] = dcd_cnts_df["truth"] - dcd_cnts_df["decoded"]

    fp = -np.sum(dcd_cnts_df["diff"].loc[dcd_cnts_df["diff"] <= 0])
    trth_cnts_df["decoded"] = dcd_cnts
    fdr = fp/np.sum(dcd_cnts)


trth_cnts_df.fillna(0, inplace=True)

trth_cnts_df["diff"] = trth_cnts_df["truth"] - trth_cnts_df["decoded"]

not_decoded_true = np.sum(trth_cnts_df["diff"].loc[trth_cnts_df["diff"] >=0])
eff = (total_trth_cnts - not_decoded_true)/total_trth_cnts





sumstats["eff_truth"] = eff
sumstats["fdr_truth"] = fdr
sumstats["fp"] = fp
sumstats["not_decoded_true"] = not_decoded_true
sumstats["total_true_cnts"] = total_trth_cnts
sumstats["nbarcodes"] = snakemake.wildcards['nbarcodes']
sumstats["rep"] = snakemake.wildcards["rep"]
sumstats["code"] = snakemake.wildcards["code"]
sumstats["ncodewords"] = len(cb) #snakemake.wildcards["ncws"]
#sumstats["nnonspec"] = snakemake.wildcards["nnonspec"]
sumstats["pdrop"] = snakemake.wildcards["pdrop"]


sumstats.to_csv(snakemake.output[0])