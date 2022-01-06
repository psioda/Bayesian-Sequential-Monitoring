#!/bin/bash

#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 07:59:59
#SBATCH --mem=1g
#SBATCH --output=../slurm/%a.slurm
#SBATCH --error=../error/%a.err
#SBATCH --array=1-800


## add R module
module add r/3.6.0

## run R command
R CMD BATCH "--no-save --args $SLURM_ARRAY_TASK_ID" 07_code_main.R ../logs/log$SLURM_ARRAY_TASK_ID.Rout
