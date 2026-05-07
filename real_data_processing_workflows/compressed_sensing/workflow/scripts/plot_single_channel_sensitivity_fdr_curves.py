import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

all_pnts = pd.concat(map(pd.read_csv, snakemake.input.best_curve_pnts))
codebook = pd.read_csv(snakemake.input.codebook)

smFISH_reference = pd.read_csv(snakemake.input.smFISH_reference)

gene_counts_by_target_fdr = all_pnts.groupby(['target_eFDR', 'gene']).size().unstack()

eFDRs = []
sensitivitites = []

for target_eFDR_results in gene_counts_by_target_fdr.groupy('target_eFDR'):
    target_eFDR, gene_counts = target_eFDR_results
    gene_counts = gene_counts.reindex(codebook.gene, fill_value=0)
    gene_counts = gene_counts.loc[smFISH_reference.gene]
    smFISH_ave_counts = smFISH_reference.ave_counts.values
    ncells_seqFISH = len(all_pnts.groupby(['pos', 'cell']))   #all_pnts.rep.nunique() * all_pnts.pos.nunique()
    
    results = pd.DataFrame({"seqFISH_ave_counts": gene_counts})/ncells_seqFISH
    
    results["smFISH_ave_counts"] = smFISH_ave_counts 
    
    results.dropna(inplace=True)
    
    lm_synd_pub = stats.linregress(results['smFISH_ave_counts'], results['seqFISH_ave_counts'])
    eFDRs.append(target_eFDR)
    sensitivitites.append(lm_synd_pub.slope)
    
    plt.scatter(results['smFISH_ave_counts'], results['seqFISH_ave_counts'], label=f"eFDR={target_eFDR:.2f}\nSensitivity={lm_synd_pub.slope:.2f}\nIntercept={lm_synd_pub.intercept:.2f}")

plt.plot(eFDRS, sensitivitites)
plt.xlabel("Estimated FDR")
plt.ylabel("Sensitivity relative to smFISH")
plt.savefig(snakemake.output[0])