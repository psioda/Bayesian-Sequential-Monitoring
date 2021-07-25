#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 15:00:00
#SBATCH --mem=1g
#SBATCH --output=../slurm/slurm.%a
#SBATCH --error=../error/%a.err
#SBATCH --array=1-100


## add R module
module add r/3.6.0

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" code_source.R ../logs/log$SLURM_ARRAY_TASK_ID.Rout
