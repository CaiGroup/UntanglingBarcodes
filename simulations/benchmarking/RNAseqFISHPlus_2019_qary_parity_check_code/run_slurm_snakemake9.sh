#!/bin/bash

#Submit this script with:sbatch thefilename

#SBATCH --time=150:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=9G   # memory per CPU core
#SBATCH -J "snakemaster"   # job name


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load julia/1.12.2
module load gurobi/10.0.0

# Snakemake 8+ moved generic cluster submission (the old --cluster flag) out
# of the core CLI and into the cluster-generic executor plugin. Install it
# once per environment with:
#   pip install snakemake-executor-plugin-cluster-generic
# --parsable makes sbatch print just the numeric job ID on stdout, which the
# plugin needs to track job status (plain sbatch's "Submitted batch job N"
# text is not parsed).
#
# The submit command templates in {resources.mem_mb}, but only one rule
# (install_cs_julia_environment) sets mem_mb explicitly. --default-resources
# fills in mem_mb for every other rule so the template always resolves.
snakemake --use-conda --cores=1 \
    --default-resources mem_mb=3000 \
    --executor cluster-generic \
    --cluster-generic-submit-cmd 'sbatch --parsable -t 5000 -c 1 --mem-per-cpu={resources.mem_mb}' \
    -j 150 --scheduler greedy --rerun-incomplete --rerun-triggers mtime --retries 3
