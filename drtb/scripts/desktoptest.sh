#!/bin/bash
# R --slave --vanilla --args < ../R/2stanfitting.R together_ARk_2 & R --slave --vanilla --args < ../R/2stanfitting.R together_ARk_4
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
rg=${key[$1-1]} # get option
echo $rg
# R --slave --vanilla --args < ../R/2stanfitting.R $rg
