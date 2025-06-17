# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 10:01:13 2021

@author: jonat
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

#res_opt = pd.read_csv("decode_summary_stats.csv")
#res_opt = pd.read_csv("decode_summary_stats_mw600.csv")

#res_opt = res_opt.loc[4:]
#res_mvc = pd.read_csv("scan_mvc_summary.csv")
#sum_stats = pd.read_csv('summary_stats_segmented_cells.csv')

sum_stats = pd.read_csv(snakemake.input[0])



rw = pd.DataFrame()
#rw["Dot Detection & Alignment"] = ["ADCG + In Channel Fiducial Markers"]*len(sum_stats['dr']) + ["LOG/Radial Center + Dapi Cross Correlation"]*(len(res_mvc["radius"]))
#rw['Decoder'] = ['Syndrome Decoding + Optimization']*(len(sum_stats['dr'])) + ["Seed Search + Minimum Variance Consensus"]*len(res_mvc["radius"])
#rw['Workflow'] = ['New']*(len(sum_stats['dr'])) + ["Published"]*len(res_mvc["radius"])

rw['Gene Encoding Barcodes (counts)'] = sum_stats['gene_encoding_barcodes']
rw['Negative Control Barcodes (counts)'] = sum_stats['negative_control_barcodes']
rw["Lateral Variance Penalty"] = sum_stats['lvf']
#rw["Number of Allowed Drops"] = sum_stats['dr']
rw["Z Variance Penalty"] = sum_stats['zvf']
rw["Log(weight) Variance Penalty"] = sum_stats['lwvf']
#rw["Estimated False Postive Rate"] = res_opt['est_fp_rate']

#sns.relplot(data=rw, x="Off Target Barcodes (counts)", y = "On Target Barcodes (counts)", hue="Decoder", style="Dot Detection & Alignment", size="Number of Allowed Drops")
sns.relplot(data=rw, x="Negative Control Barcodes (counts)", y = "Gene Encoding Barcodes (counts)", hue="Lateral Variance Penalty", size="Z Variance Penalty", style = "Log(weight) Variance Penalty")

plt.tight_layout()
plt.savefig(snakemake.output[0])
#plt.show()
#sns.relplot(data=res, x="Estimated False Postive Rate", y = "On Target Barcodes (counts)", hue="Lateral Variance Penalty", size="Number of Allowed Drops")


#rw = rw[rw['dr'] == 0]

#rw["log2(weight) Varience Penalty"] = rw['lwvf']
#rw['Sigma Varience Penalty'] = rw['svf']

#sns.relplot(data=rw, x="Estimated False Postive Rate", y = "On Target Barcodes (counts)", hue="Lateral Variance Penalty", style="log2(weight) Varience Penalty", size = 'Sigma Varience Penalty')

#rw = rw[rw['svf'] == 0]

#sns.relplot(data=rw, x="Off Target Barcodes (counts)", y = "On Target Barcodes (counts)", hue="Lateral Variance Penalty", style="log2-weight Varience Penalty")
#sns.relplot(data=rw, x="Estimated False Postive Rate", y = "On Target Barcodes (counts)", hue="Lateral Variance Penalty", style="log2(weight) Varience Penalty")
