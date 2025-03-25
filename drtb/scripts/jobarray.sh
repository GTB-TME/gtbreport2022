#!/bin/bash
#SBATCH --job-name=STOUT
#SBATCH --time=7:59:00
#SBATCH --mem=10G
#SBATCH --output=test_%A_%a.out
#SBATCH --error=test_%A_%a.err
#SBATCH --array=1-15
module load R/4.0.0-foss-2020a

# # root names:
# separate_LerouxIntercept
# together_ARk
# separate_LerouxInterceptLerouxSlope
# together_LerouxIntercept
# separate_global
# together_LerouxInterceptLerouxSlope
# separate_globalH

key=(separate_LerouxIntercept_2
     together_ARk_2
     separate_LerouxInterceptLerouxSlope_2
     together_LerouxIntercept_2
     separate_global_2
     together_LerouxInterceptLerouxSlope_2
     separate_globalH_2
     # 4 versions
     separate_LerouxIntercept_4
     together_ARk_4
     separate_LerouxInterceptLerouxSlope_4
     together_LerouxIntercept_4
     separate_global_4
     together_LerouxInterceptLerouxSlope_4
     separate_globalH_4
) # options key NOTE change array length!
rg=${key[$SLURM_ARRAY_TASK_ID-1]} # get option
R --slave --vanilla --args < ../R/2stanfitting.R $rg

