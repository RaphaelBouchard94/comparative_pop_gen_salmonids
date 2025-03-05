rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)


#---------------#
#----- Data ----#
#---------------#

#Whitefish

cocl_zscore <- read.table("data/wgs_assign/whitefish/cocl_zscore_site.txt",header=F)

colnames(cocl_zscore) <- c("id","site","z")

u_var_per_site <- cocl_zscore %>% group_by(site) %>% 
  summarise(mean_z = mean(z),
            var_z = var(z))

cocl_zscore %<>% left_join(u_var_per_site, by = "site")

cocl_zscore %<>% mutate(lk = (z-mean_z)/var_z)

cocl_zscore %>% group_by(site) %>% mutate(sd_lk = sd(lk))

pop_cocl <- unique(cocl_zcore$site)

zdeviant_cocl <- list()

for(i in 1:length(pop)) {
tresh <- quantile(subset(cocl_zscore,subset = site == pop_cocl[i])$lk,probs = c(0.1,0.9))
mig <- cocl_zscore %>% filter(site == pop_cocl[i], lk > tresh[2] | lk < tresh[1]) %>% pull(id)
zdeviant_cocl[[i]] <- mig
}

migrant_id_cocl <- unlist(zdeviant_cocl)

cocl_zscore %>% ggplot(aes(x = site, y = lk, color = id %in% migrant_id_cocl))+
  geom_jitter(width = 0.1,alpha = 0.6)+
  scale_color_manual(values = c("black","red"))+
  theme_bw()


#Trout

safo_zscore <- read.table("data/wgs_assign/trout/safo_zscore.txt",header=F)

colnames(safo_zscore) <- c("id","site","z")

u_var_per_site <- safo_zscore %>% group_by(site) %>% 
  summarise(mean_z = mean(z),
            var_z = var(z))

safo_zscore %<>% left_join(u_var_per_site, by = "site")

safo_zscore %<>% mutate(lk = (z-mean_z)/var_z)


safo_zscore %>% ggplot(aes(x = site, y = lk))+
  geom_jitter(width = 0.1,alpha = 0.6)

pop <- unique(safo_zscore$site)

zdeviant <- list()

for(i in 1:length(pop)) {
  tresh <- quantile(subset(safo_zscore,subset = site == pop[i])$lk,probs = c(0.1,0.9))
  mig <- safo_zscore %>% filter(site == pop[i], lk > tresh[2] | lk < tresh[1]) %>% pull(id)
  zdeviant[[i]] <- mig
}

migrant_id <- unlist(zdeviant)

safo_zscore %>% ggplot(aes(x = site, y = lk, color = id %in% migrant_id))+
  geom_jitter(width = 0.1,alpha = 0.6)+
  scale_color_manual(values = c("black","red"))+
  theme_bw()



