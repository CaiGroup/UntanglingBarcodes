# Codebook generation

To generate a codebook for a seqFISH experiment, you first must generate a table of codewords from an error correcting codes. More information on this step can be found in the [get_RS_codebooks](https://github.com/CaiGroup/UntanglingBarcodes/tree/main/codebook_generation/get_RS_codebooks) readme and its accompanying notebooks.

## Installing a julia kernal to run these notebooks in Jupyter
To run these notebooks, you will need to install a Julia kernel. If you are running the notebooks localling, first install Julia by downloading it from the julia [website](https://julialang.org). If you are on the Caltech HPC, load a pre-installed julia module by typing

```
module load julia/1.10.8
```
into your bash terminal.

Next open a Julia session in a terminal and install the IJulia package by typing
```
julia> using Pkg
julia> Pkg.add("IJulia")
julia> Pkg.build("IJulia")
```

After the package installation is complete, a Julia kernel should be available in Jupyter.
