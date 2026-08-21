import pandas as pd

dfs = [pd.read_csv(f) for f in snakemake.input]
result = pd.concat(dfs, ignore_index=True)
result.to_csv(snakemake.output[0], index=False)
