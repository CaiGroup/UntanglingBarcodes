# Relative code robustness simulation

To run this simulation, first you will need to install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) 7, the version of snakemake in which this workflow was tested. I cannot guarantee that it will work in other versions of snakemake.

To install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) 7, use the following commands: 

```
conda activate base
conda create -c conda-forge -c bioconda -n snakemake7 snakemake=7
conda activate snakemake7
```

You will need to activate your snakemake conda environment whenever you run the workflow.

Before running the workflow, it is necessary to install the Julia pacakges:

DataFrames.jl
CSV.jl
Gurobi.jl
SeqFISHSyndromeDecoding.jl
Distributions.jl

You can install these packages by opening a julia session and typing the commands

```
>>using Pkg
>>Pkg.add("CSV")
>>Pkg.add("DataFrames")
>>Pkg.add(name="Gurobi",version="1.2.3")
>>Pkg.add("Distributions")
```

To install SeqFISHSyndromeDecoding.jl, navigate into the package directory with the julia session, then run the command

```
>>using Pkg
>>Pkg.add(".")
```

Snakemake will automatically install all necessary python packages on its own on linux or mac. In windows run the following commands on command prompt install the following packages into your snakemake conda environment:

```
conda install numpy
conda install pandas
```

Once the dependencies are installed, run the simulation locally with

```
snakemake --use-conda -c<n>
```
where <n> is the number of cores you wish to provide. On a slurm cluster, submit the batch script using the command

```
sbatch run_slurm.sh
```