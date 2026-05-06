import pandas as pd

dfs = []
for file in snakemake.input:
    hybdf = pd.read_csv(file)
    hybdf['hyb'] = int(file.split('.')[0].split('hyb_')[1])
    dfs.append(hybdf)

df = pd.concat(dfs)
df.set_index('hyb', inplace=True)
df.to_csv(snakemake.output[0])
