rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#-----------------#
#---- LIBRARY ----#
#-----------------#

library(tidyverse)
library(magrittr)
library(readxl)
library(patchwork)

#-----------------#
#------ DATA -----#
#-----------------#

safo_ne <- read_xlsx("data/ne/safo_ldne.xlsx")

cocl_ne <- read_xlsx("data/ne/cocl_ldne.xlsx")

ne <- rbind(safo_ne,cocl_ne)

ne$sp <- c(rep("S. fontinalis",10),rep("C. clupeaformis",8))

sites <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/pop_dispersal.txt",header=T)

ne <- left_join(ne,sites, by = c("pop"="pop_code"))

ne %<>% mutate(pop2 = case_when(pop == "KAA" & sp == "S. fontinalis" ~ "ROG-K",
                               pop == "LAG" & sp == "C. clupeaformis" ~ "BIA-L",
                               pop == "EAS" & sp == "C. clupeaformis" ~ "EAS-C",
                               T ~ pop))

#-----------------#
#---- ANALYSIS ---#
#-----------------#

ne %>% ggplot(aes(x = ne,y = reorder(pop2,lat), shape = sp))+
  geom_point(position = position_dodge(.5), size= 1)+
  geom_errorbar(aes(xmin=lower_ci,xmax=upper_ci,linetype = sp),width = 0,position = position_dodge(.5), size = .3)+
  scale_linetype_manual(values = c(2,1), labels = c(expression(italic("C. clupeaformis")),expression(italic("S. fontinalis"))))+
  scale_shape_manual(values = c(1,19),labels = c(expression(italic("C. clupeaformis")),expression(italic("S. fontinalis"))))+
  ylab("Population")+
  xlab("Estimated Ne value")+
  coord_cartesian(xlim = c(0, 2000)) +  
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18,hjust = 0),
        legend.position = c(0.8,0.9),
        legend.background = element_rect(color = "black"))

ggsave("~/Desktop/figures_paper/ne_per_pop.png")

ne %>% group_by(sp) %>% summarize(median_ne = median(ne))

t.test(x = ne %>% filter(sp == "C. clupeaformis") %>% dplyr::select(ne), y = ne %>% filter(sp == "S. fontinalis") %>% dplyr::select(ne))

ne %>% ggplot(aes(x = ne,y = reorder(pop2,lat)))+
  geom_point(position = position_dodge(.5), size= 1.5)+
  geom_errorbar(aes(xmin=lower_ci,xmax=upper_ci),width = 0,position = position_dodge(.5), size = .3)+
  scale_linetype_manual(values = c(2,1))+
  scale_shape_manual(values = c(1,19))+
  ylab("Population")+
  xlab("Estimated Ne value")+
  coord_cartesian(xlim = c(0, 2000))+
  facet_wrap(~ sp, nrow = 2, ncol = 1,scales = "free_y")+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18,hjust = 0))

ne_cocl <- ne %>% filter(sp == "C. clupeaformis") %>% ggplot(aes(x = ne,y = reorder(pop2,lat)))+
  geom_point(position = position_dodge(.5), size= 1)+
  geom_errorbar(aes(xmin=lower_ci,xmax=upper_ci),width = 0.2,position = position_dodge(.5), size = .3)+
  ylab("Population")+
  xlab("")+
  coord_cartesian(xlim = c(0,2000))+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,hjust = .8),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18,hjust = 0))

ne_cocl

ne_safo <- ne %>% filter(sp == "S. fontinalis") %>% ggplot(aes(x = ne,y = reorder(pop2,lat)))+
  geom_point(position = position_dodge(.5), size= 1)+
  geom_errorbar(aes(xmin=lower_ci,xmax=upper_ci),width = 0.2,position = position_dodge(.5), size = .3)+
  ylab("Population")+
  xlab("Estimated Ne value")+
  coord_cartesian(xlim = c(0,2000))+
  theme_bw()+
  theme(strip.background = element_rect(fill=NA),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18,hjust = 0))


ne_cocl + whitefish2 + ne_safo + trout2

#ggsave("~/Desktop/figures_paper/ne_het_per_pop.png",width = 9, height = 8)

lat_ne <- lm(ne ~ sp + lat, data = ne)

summary(lat_ne)

plot(lat_ne)

#
ne %>% group_by(sp) %>% summarise(median_n = median(ne),
                                  mean_n = mean(ne))




