import pandas as pd

reps = [pd.read_csv(f) for f in snakemake.input]

reps_df = pd.concat(reps)

reps_grp_sum = reps_df.groupby(['lvf','lwvf','svf','dr']).sum()

reps_grp_sum.drop(['pos','ch','est_fp_rate'], axis="columns")


n_gene_encoding_cws = 10000
n_neg_ctrl_cws = 14000

reps_grp_sum['est_fp_rate'] = n_gene_encoding_cws*(reps_grp_sum['negative_control_barcodes']/n_neg_ctrl_cws)/reps_grp_sum['gene_encoding_barcodes']

reps_grp_sum.to_csv(snakemake.output[0])
