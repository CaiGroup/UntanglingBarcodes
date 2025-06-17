import pandas as pd

barcodes = pd.read_csv(snakemake.input[0])
replicate = int(snakemake.wildcards["rep"])

for infile in snakemake.input[1:]:
    barcodes = barcodes.append(pd.read_csv(infile))

barcodes = barcodes.groupby(["lvf", "wf", "svf", "dr", "pos", "ch", "cellid", "gene"]).size()

#barcodes.rename(columns={"0":"counts"}, inplace=True)
barcodes.name="counts"
barcodes["replicate"] = replicate
barcodes.to_csv(snakemake.output[0], index=True)
