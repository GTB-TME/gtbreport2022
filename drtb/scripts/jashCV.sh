#!/bin/bash
#$ -N CARfitCV4
#$ -cwd -V
#$ -l h_rt=7:59:00
#$ -l rmem=10G
#$ -t 1-97
module load apps/R/4.0.0/gcc-10.1
R --slave --vanilla --args < ../R/2stanfitting.R together_LerouxInterceptLerouxSlope_2 $SGE_TASK_ID
# R --slave --vanilla --args < ../R/2stanfitting.R separate_LerouxInterceptLerouxSlope_2 $SGE_TASK_ID
# R --slave --vanilla --args < ../R/2stanfitting.R together_LerouxIntercept_4 $SGE_TASK_ID
# R --slave --vanilla --args < ../R/2stanfitting.R separate_global_4 $SGE_TASK_ID

