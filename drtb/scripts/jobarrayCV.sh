#!/bin/bash
#SBATCH --job-name=STOUT-CV
#SBATCH --time=7:59:00
#SBATCH --mem=10G
#SBATCH --output=test_%A_%a.out
#SBATCH --error=test_%A_%a.err
#SBATCH --array=1-39
module load R/4.0.0-foss-2020a
R --slave --vanilla --args < ../R/2stanfitting.R together_LerouxInterceptLerouxSlope_2 $SLURM_ARRAY_TASK_ID
