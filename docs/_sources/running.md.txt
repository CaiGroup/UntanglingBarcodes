# Running

In this page, I give instructions on how to run the workflow locally or on a slurm cluster. The default configuration is to run all fields of view. If you want 
to test the workflow on a small portion of the data, edit the positions lists in the config file.

For the Reed-Solomon encoded processing workflow, edit the python list in line 23

```positions: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30] ```

to include only the positions that you would like to test on.

For the seqFISH+ 2019 experiment workflow, edit the replicate position lists in lines 26 and 27

```
positions:
  rep1: [0,1,2,3,4,5,6]
  rep2: [0,1,2,3,4,5,6,7,8,9]
```
or the channel list in line 19

```channels: [488, 561, 643]```

The workflow requires that at least one position be processed in each replicate, or it will error. The channels are decoded independently, so not all need to be decoded.

## Local

To run locally open your conda environment with snakemake installed in it (from above), navigate to the repository home directory, then use the command:

<pre> <code> snakemake --use-conda --cores n </code> </pre>

where n is the number of cores you would like to provide.

## Slurm Cluster (as in the Caltech HPC)

To run on the HPC open your conda environment with snakemake installed in it (from above), then run the command:

<pre> <code> sbatch run_slurm.sh </code> </pre>

from the repository home directory.

If you have a large dataset, the master job may time out. To increase the allowed time, change edit line 5 of run-slurm.sh:

<pre> <code> #SBATCH --time=150:00:00   # walltime </code> </pre>

I have set the walltime to 150 hours, which I find is enough for my datasets, but the Caltech HPC allows walltimes of up to 14 days.

You may also want to change the maximum number of jobs that may be submitted at a time. To change this, edit line 16 of run-slurm.sh:

<pre> <code> snakemake --use-conda --cores=1 --cluster 'sbatch -t 3000 -c 1' -j 150 --scheduler greedy --rerun-incomplete --keep-going --rerun-tiggers mtime</code> </pre>

The last option, <code> -j </code>, specifies the maximum number of jobs that may be submitted at a time. I have it set to 150.

I have found that error occasionally occur with submitting many parallel jobs or with requesting Gurobi licenses from multiple parallel jobs. These errors may prevent
the workflow from finishing in one run. If they do occur, just restart the workflow after each error and all jobs should eventually complete successfully.

## Useful Commands

To check how many jobs are left to execute, run

<pre> <code> snakemake -n </code> </pre>

Some useful options for running worfkows with the 'snakemake' command

- `-scheduler greedy` prevents snakemake from spending too much time optimizing the order in which to submit jobs
- `--keep-going` instructs it to evaluate all jobs whose input files can be generated even if there are errors in other jobs.
- `--rerun-tiggers mtime` is useful to avoid needlessly rerunning jobs.
- `--use-conda` instructs snakemake to install conda environments for each job (does not work on windows)

### Useful Slurm Commands

To see what slurm jobs you have running on the HPC use the command

<pre> <code> squeue -u your_username </code> </pre>

To check for errors, search the output files with

<pre> <code> grep -i error ./* </code> </pre>
