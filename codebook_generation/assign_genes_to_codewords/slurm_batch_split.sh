#!/bin/bash

#Submit this script with:sbatch thefilename

#SBATCH --time=150:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=10G   # memory per CPU core
#SBATCH -J "split codebooks"   # job name


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load julia/1.10.8
module load gurobi/10.0.0

julia codebook_split.jl