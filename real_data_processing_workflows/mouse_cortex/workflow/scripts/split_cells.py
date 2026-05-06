import pandas as pd

fov_dots = pd.read_csv(snakemake.input[0])


for outfname in snakemake.output:
    print(outfname)
    cell_dot_diff = outfname.split("_cell")[1]
    cell = int(cell_dot_diff.split(".")[0])
    print("cell: ", cell)
    cell_dots = fov_dots.loc[fov_dots.cellid == cell]
    cell_dots.to_csv(outfname)
