rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(plyr)
library(RcppCNPy)


#---------------#
#----- Data ----#
#---------------#


af_trout <- npyLoad("data/wgs_assign/trout/safo_254_essbp_ind_118861_pruned_singleton_snp.pop_af.npy")
