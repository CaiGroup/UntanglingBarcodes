# RS SeqFISH workflow
This is the workflow used to process raw data from our Reed-Solomon encoded experiment. The raw images are too large to include in this repository, but an alternative version of the workflow is available in another [git repository](https://github.com/CaiGroup/DisentanglingRSBarcodes/tree/master) that includes preprocessed data of a managable size. 


# Installing packages

To run this simulation, first you will need to install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) 7, the version of snakemake in which this workflow was tested. I cannot guarantee that it will work in other versions of snakemake.

To install [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) 7, use the following commands: 

```
conda activate base
conda create -c conda-forge -c bioconda -n snakemake7 snakemake=7
conda activate snakemake7
```

You will need to activate your snakemake7 conda environment whenever you run the workflow.

Before running the workflow, it is necessary to install the Julia pacakges:

DataFrames.jl
CSV.jl
Images.jl
FileIO.jl
Gurobi.jl
SeqFISH_ADCG.jl
SeqFISHSyndromeDecoding.jl

You can install these packages by opening a julia session and typing the commands

```
>>using Pkg
>>Pkg.add("CSV")
>>Pkg.add("DataFrames")
>>Pkg.add("Image")
>>Pkg.add("FileIO")
>>Pkg.add("https://github.com/CaiGroup/SeqFISH_ADCG.jl")
>>Pkg.add("https://github.com/CaiGroup/SeqFISHSyndromeDecoding.jl")
>>Pkg.add(name="Gurobi",version="1.2.3")
```

Snakemake will automatically install all necessary python packages on its own on mac or linux operating systems.


To run the workflow locally, use the command
```
snakemake --use-conda -c<n>
```
where `<n>` is the number of cores that you would like to assign to the workflow.

