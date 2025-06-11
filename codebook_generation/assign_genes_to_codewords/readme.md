# Assigning genes to codewords

After choosing to use weight 4 codewords of the q11n10k7 Reed-Solomon code and finding the codewords, I designed the experimental implementation in a few steps. First I chose a 542 codewords subset of the q11n10k7 Reed-Solomon weight 4 codewords using `codebook_thinning.jl` that maximizes the Hamming distances between codewords. Then I assigned genes to the codewords using `codebook_equalization.jl`. This script takes as input the table of codewords output from `codebook_thinning.jl` and a table of genes with their bulk RNAseq expression levels, then assigns genes to codewords to minimize the maximum sum of RNAseq expression levels all genes probed in any single hybridization cycle. Finally, I assigned block-pseudocolors to hybridization-channels using `intra_channel_optimization.jl` and split the experiment into two lower density pools using `codebook_split.jl`.

I do not run these optimizations to find the absolute best solution, but end them after either reaching a time limit or finding a solution within an acceptable margin of the bound on the optimal solution, so I cannot guarantee that you will get the same results from these scripts as I did. For demonstration purposes, I have included intermediate files I used in designing the experiments in the resources in the resources folder.

The `Project.toml` and `Manifest.toml` files include the julia environment to run these scripts.

To run these optimizations, I use [JuMP](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) modeling language with the [Gurobi](https://www.gurobi.com/) commercial integer linear programming solver. Alternative [open source solvers](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers) are available, but are slower. To use these scripts you will need to either install Gurobi or modify the scripts to use an open source solver.

## Running scripts Locally

Once Gurobi is installed, or the scripts have been edited to use another solver which has been installed, the scripts can be run as normal julia scripts, by entering

```julia <script name>.jl```

## Running scripts on a slurm cluster

I include slurm batch submission scripts (the files with names of the form slurm_batch_*.sh) for running each of the julia scripts on a slurm cluster. 

For example, to submit the a batch job running the thin_codebooks script on a slurm cluster, navigate to this folder in a shell window and then type

``` sbatch slurm_batch_thin.sh ```

The Caltech HPC has gurobi installed as a module and can be loaded by adding the line. If you use a different slurm cluster, you may need to modify the Gurobi version loaded in the batch script or use an open source solver.
