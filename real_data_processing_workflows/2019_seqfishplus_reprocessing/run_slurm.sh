#!/bin/bash

#Submit this script with:sbatch thefilename

#SBATCH --time=124:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=4G   # memory per CPU core
#SBATCH -J "snakemaster"   # job name


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load julia/1.10.2
module load gurobi/9.5.1

snakemake --use-conda --cores=1 --cluster 'sbatch -t 5000 -c 1 --mem-per-cpu={resources.mem_mb}' -j 30 --scheduler greedy --rerun-incomplete --keep-going --rerun-triggers mtime
