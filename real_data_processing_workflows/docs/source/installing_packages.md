# Installing packages

To install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html), use the following commands from the snakemake full installation [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html): 

```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```


The [SeqFISHSyndromeDecoding.jl](https://github.com/CaiGroup/SeqFISHSyndromeDecoding) package uses an integer programming solver. In our run of the workflow, we used a commercial solver, [Gurobi](https://www.gurobi.com/), which provides free academic licenses.
Using the workflows as written requires that Gurobi be installed. Alternatively, the `choose_cpaths_from_saved_candidates.jl` and `choose_cpaths_from_saved_candidates_no_neg_ctrl.jl` scripts may be modified to use another [integer linear programming solver](https://jump.dev/JuMP.jl/stable/installation/#Supported-solvers). Each of these scripts contains a commented alternative funtion call to `SeqFISHSyndromeDecoding.choose_optimal_codepaths` using the opensource GLPK integer linear programming solver.


You will need to activate your snakemake conda environment whenever you run the workflow.

Before running the workflow, it is necessary to install the Julia pacakges:

DataFrames.jl
CSV.jl
Images.jl
FileIO.jl
JuMP.jl
Gurobi.jl
SeqFISH_ADCG.jl
SeqFISHSyndromeDecoding.jl

You can install these packages by opening a julia session and typing the commands

```
>>using Pkg
>>Pkg.add("CSV")
>>Pkg.add("DataFrames")
>>Pkg.add("Images")
>>Pkg.add("FileIO")
>>Pkg.add("JuMP")
>>Pkg.add("Gurobi")
>>Pkg.add("https://github.com/CaiGroup/SeqFISH_ADCG")
>>Pkg.add("https://github.com/CaiGroup/SeqFISHSyndromeDecoding")
```

Snakemake will automatically install all necessary python packages on its own in linux or max operating systems. The python packages listed in each .yaml file of the workflow/envs folder will need to be installed manually with conda to run the workflows on windows machines.