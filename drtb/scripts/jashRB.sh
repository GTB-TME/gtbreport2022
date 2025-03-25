#!/bin/bash
#$ -N CARfitRB
#$ -cwd -V
#$ -l h_rt=7:59:00
#$ -l rmem=10G
#$ -t 1-5
module load apps/R/4.0.0/gcc-10.1
R --slave --vanilla --args < ../R/2stanfitting.R together_LerouxInterceptLerouxSlope_2 $((1-$SGE_TASK_ID)) T

