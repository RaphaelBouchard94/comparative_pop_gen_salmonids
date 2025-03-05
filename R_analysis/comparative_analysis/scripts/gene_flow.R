rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(marmap)
library(reshape2)

#---------------#
#----- Data ----#
#---------------#

####################################
######## Bathymetric map  ##########
####################################

JB <- getNOAA.bathy(lon1 = -82.084, lon2 = -78,
                    lat1 = 51.021, lat2 = 55, resolution = 1)

sites <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/pop_dispersal.txt",header=T)

marsites<- sites %>% select(long,lat) %>% rename(x=long,y=lat)

trans1 <- trans.mat(JB,min.depth = 5)

dist1 <- lc.dist(trans1,marsites,res="dist")

dist1

mat_dist<-as.matrix(dist1)

dfdist <- setNames(melt(mat_dist), c('pop1', 'pop2', 'values'))

###################
###  Gene flow ####
###################

fst_cocl <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/fst/fst_cocl.tsv")

colnames(fst_cocl) <- c("fst","pop1","pop2","m")

fst_cocl$species <- rep("Coregonus clupeaformis",78)

fst_safo <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/fst/fst_safo.tsv")

colnames(fst_safo) <- c("fst","pop1","pop2","m")

fst_safo$species <- rep("Salvelinus fontinalis",91)

fst <- rbind(fst_cocl,fst_safo)

fst$pop1 <- gsub("CHI","KAA",fst$pop1)
fst$pop2 <- gsub("CHI","KAA",fst$pop2)
fst$pop1 <- gsub("FIS","EAS",fst$pop1)
fst$pop2 <- gsub("FIS","EAS",fst$pop2)
fst$pop1 <- gsub("MOA","SAB",fst$pop1)
fst$pop2 <- gsub("MOA","SAB",fst$pop2)
fst$pop1 <- gsub("WAS","JAC",fst$pop1)
fst$pop2 <- gsub("WAS","JAC",fst$pop2)

fst$pop1 <- factor(fst$pop1, levels=c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","JAC","RUP","BRO","NOT"))
fst$pop2 <- factor(fst$pop2, levels=c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","JAC","RUP","BRO","NOT"))


fst %>% ggplot(aes(x = species, y = log10(m)))+
  geom_boxplot(width = 0.2)+
  geom_jitter(alpha = 0.4,width = 0.03)+
  ylab("log10(Gene flow)")+
  xlab("")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))



dfdist %<>% mutate(pop1_code = case_when(pop1 == 1 ~ "ROG",
                                         pop1 == 2 ~ "KAA",
                                         pop1 == 3 ~ "BIA",
                                         pop1 == 4 ~ "LAG",
                                         pop1 == 5 ~ "AQU",
                                         pop1 == 6 ~ "BEA",
                                         pop1 == 7 ~ "SAB",
                                         pop1 == 8 ~ "REN",
                                         pop1 == 9 ~ "CON",
                                         pop1 == 10 ~ "EAS",
                                         pop1 == 11 ~ "JAC",
                                         pop1 == 12 ~ "RUP",
                                         pop1 == 13 ~ "BRO",
                                         pop1 == 14 ~ "NOT"),
                   pop2_code = case_when(pop2 == 1 ~ "ROG",
                                         pop2 == 2 ~ "KAA",
                                         pop2 == 3 ~ "BIA",
                                         pop2 == 4 ~ "LAG",
                                         pop2 == 5 ~ "AQU",
                                         pop2 == 6 ~ "BEA",
                                         pop2 == 7 ~ "SAB",
                                         pop2 == 8 ~ "REN",
                                         pop2 == 9 ~ "CON",
                                         pop2 == 10 ~ "EAS",
                                         pop2 == 11 ~ "JAC",
                                         pop2 == 12 ~ "RUP",
                                         pop2 == 13 ~ "BRO",
                                         pop2 == 14 ~ "NOT")) %>% select(-(pop1:pop2))

fst %<>% left_join(dfdist,by=c("pop1"="pop1_code","pop2"="pop2_code"),relationship = "many-to-many")


fst %<>% mutate(label_comp = paste(pop1,pop2,sep="_"))

ggplot(fst,aes(x = values, y = log10(m), shape=species, linetype = species))+
  geom_point(size=3,alpha = 0.8)+
  geom_smooth(method="lm",color="black",alpha = 0.3,size = 1,se = 0.95)+
  scale_shape_manual(values = c(16,1),labels = c(expression(italic("C. clupeaformis")),expression(italic("S. fontinalis"))))+
  scale_linetype_manual(values = c(2,3),labels = c(expression(italic("C. clupeaformis")),expression(italic("S. fontinalis"))))+
  xlab("Geographic distance (km)")+
  ylab("log10(Gene flow)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = c(0.8,0.9),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))

ggsave("figures/gene_flow.jpeg",width = 8, height = 8)



fst %>% group_by(species) %>% summarise(avg_m = mean(m),
                                        ci_m = sd(m))

