# Installing packages

To install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html), use the following commands from the snakemake full installation [instructions](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html): 

```
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

You will need to activate your snakemake conda environment whenever you run the workflow.

Before running the workflow, it is necessary to install the Julia pacakges:

DataFrames.jl
CSV.jl
Images.jl
FileIO.jl
SeqFISH_ADCG.jl
SeqFISHSyndromeDecoding.jl

You can install these packages by opening a julia session and typing the commands

```
>>using Pkg
>>Pkg.add("CSV")
>>Pkg.add("DataFrames")
>>Pkg.add("Image")
>>Pkg.add("FileIO")
>>Pkg.add("https://github.com/CaiGroup/SeqFISH_ADCG")
>>Pkg.add("https://github.com/CaiGroup/SeqFISHSyndromeDecoding")
```

Snakemake will automatically install all necessary python packages on its own.