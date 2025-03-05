rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(devtools)
source("scripts/plotting_funcs.R")

#---------------#
#----- Data ----#
#---------------#

plot_tree("data/treemix/cocl/m/cocl.treemix.m1.rep1", )

plot_resid("data/treemix/cocl/m/cocl.treemix.m1.rep1",pop_order = "data/treemix/cocl/pop_order.txt")



pdf("figures/treemix_coclm1.pdf",width = 7,height = 7)
plot_tree("data/treemix/cocl/m/cocl.treemix.m1.rep9",mbar =F,scale = F)
dev.off()

plot_tree("data/treemix/cocl/m/cocl.treemix.m0.rep10",mbar =F,scale = F)

pdf("figures/treemix_safom1.pdf",width = 7,height = 7)
plot_tree("data/treemix/safo/m/safo.treemix.m1.rep6",ybar = 0.5,flip = "yes",mbar =F)
dev.off()

pdf("figures/treemix_safom1_zoom.pdf",width = 7,height = 7)
plot_tree("data/treemix/safo/m/safo.treemix.m1.rep6",plotmig = T,xmin = 0.15)
dev.off()

plot_tree("data/treemix/safo/m/safo.treemix.m0.rep1",disp = 0.005)


