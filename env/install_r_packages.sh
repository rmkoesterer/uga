#!/bin/bash

## note: run 'source activate uga' before installing R packages
## run this from the command line to get the current package versions, otherwise same versions can be installed manually

## versions installed during development
# kinship2_1.8.4
# geepack_1.2-1
# lme4_1.1-21
# lmerTest_3.1-0
# pbkrtest_0.4-7
# seqMeta_1.6.7
# RColorBrewer_1.1-2
# R.utils_2.9.0
# ggplot2_3.2.0

R -e 'install.packages(c("kinship2", "geepack", "lme4", "lmerTest", "pbkrtest", "seqMeta", "RColorBrewer", "ggplot2","R.utils"), repos="http://cran.us.r-project.org", dependencies=TRUE)'
