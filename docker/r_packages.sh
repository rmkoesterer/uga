#!/bin/bash

conda activate uga

R -e 'install.packages(c("kinship2", "geepack", "lme4", "lmertest", "pbkrtest", "seqMeta", "RColorBrewer", "ggplot2","R.utils"), repos="http://cran.us.r-project.org", dependencies=TRUE)'
