#!/bin/bash
#$ -N CARfit
#$ -cwd -V
#$ -l h_rt=7:59:00
#$ -l rmem=10G
#$ -t 1-12
module load apps/R/4.0.0/gcc-10.1

# # root names:
# separate_LerouxIntercept
# together_ARk
# separate_LerouxInterceptLerouxSlope
# together_LerouxIntercept
# separate_global
# together_LerouxInterceptLerouxSlope
# separate_globalH

key=(separate_LerouxIntercept_2
     separate_LerouxInterceptLerouxSlope_2
     together_LerouxIntercept_2
     separate_global_2
     together_LerouxInterceptLerouxSlope_2
     separate_globalH_2
     # 4 versions
     separate_LerouxIntercept_4
     separate_LerouxInterceptLerouxSlope_4
     together_LerouxIntercept_4
     separate_global_4
     together_LerouxInterceptLerouxSlope_4
     separate_globalH_4
    ) # options key NOTE change array length!
rg=${key[$SGE_TASK_ID-1]} # get option
R --slave --vanilla --args < ../R/2stanfitting.R $rg
