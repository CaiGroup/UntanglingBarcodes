# Running

## Local

To run locally open your conda environment with snakemake installed in it (from above), navigate to the repository home directory, then use the command:

<pre> <code> snakemake --use-conda --cores n </code> </pre>

where n is the number of cores you would like to provide.

## Slurm Cluster (as in the Caltech HPC)

To run on the HPC open your conda environment with snakemake installed in it (from above), then run the command:

<pre> <code> sbatch run_slurm.sh </code> </pre>

from the repository home directory.

If you have a large dataset, the master job may time out. To increase the allowed time, change edit line 5 of run-slurm.sh:

<pre> <code> #SBATCH --time=24:00:00   # walltime </code> </pre>

The default wall time is 24 hours, but the Caltech HPC allows walltimes of up to 14 days.

You may also want to change the maximum number of jobs that may be submitted at a time. To change this, edit line 16 of run-slurm.sh:

<pre> <code> snakemake --use-conda --cores=1 --cluster 'sbatch -t 3000 -c 1' -j 200 </code> </pre>

The last option, <code> -j </code>, specifies the maximum number of jobs that may be submitted at a time. I have it set to 200. I have also tried 500, but ran in to issues with the master job not being able to manage that many jobs. The best number to use may depend on the dataset.

## Useful Commands

To check how many jobs are left to execute, run

<pre> <code> snakemake -n </code> </pre> 


### Useful Slurm Commands

To see what slurm jobs you have running on the HPC use the command

<pre> <code> squeue -u your_username </code> </pre>

To check for errors, search the output files with

<pre> <code> grep -i error ./* </code> </pre>
