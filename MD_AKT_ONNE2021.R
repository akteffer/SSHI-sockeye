## Mass Deviation analysis
# emulating Art's code in STAN

# SSHI ANALYSIS SOCKEYE STAGE II
# A.K. Teffer
#### Sockeye salmon productivity versus infection profiles

#### Load packages and set directory
#setwd("~/Documents.nosync/DFO PDF/Data/SSHI-sockeye")

library(lme4)
library(rstanarm) # https://mc-stan.org/users/documentation/case-studies/tutorial_rstanarm.html
library(ggplot2)
library(plotrix)
library(tidyverse)
library(gridExtra)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(shinystan)
library(data.table)
library(base)
library(ggpubr)
library(dplyr)

## Run code in Raw data filtering to produce final all.data df
head(all.data)


