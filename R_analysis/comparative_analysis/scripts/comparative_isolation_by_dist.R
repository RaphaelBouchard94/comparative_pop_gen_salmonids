rm(list = ls())

####################################
###########Library##################
####################################


library(tidyverse)
library(magrittr)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(readxl)
library(marmap)
library(reshape2)
library(vegan)
library(ggrepel)
library(paletteer)
library(adespatial)


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

####################################
###  FST Isolation by distance  ####
####################################

fst_cocl <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/fst/fst_cocl.tsv")

colnames(fst_cocl) <- c("fst","pop1","pop2")

fst_cocl$species <- rep("Coregonus clupeaformis",78)

fst_safo <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/fst/fst_safo.tsv")

colnames(fst_safo) <- c("fst","pop1","pop2")

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

mat_dist<-as.matrix(dist1)

dfdist <- setNames(melt(mat_dist), c('pop1', 'pop2', 'values'))

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

#Average and sd fst

fst %>% group_by(species) %>% summarise(mean_fst = mean(fst),
                                        sd_fst = sd(fst))

#Isolation by distance graph

fst %<>% mutate(label_comp = paste(pop1,pop2,sep="_"))

ggplot(fst,aes(x = values, y = fst/(1-fst)))+
  geom_point(size=3,alpha = 0.8)+
  geom_smooth(method="lm",color="black",alpha = 0.3,size = 1,linetype = 2,se = 0.95)+
  scale_color_manual(values = c("grey80","#E6452EFF"),labels = c(expression(italic("Corenus clupeaformis")),expression(italic("Salvelinus fontinalis"))))+
  scale_fill_manual(values = c("grey80","#E6452EFF"),labels = c(expression(italic("Corenus clupeaformis")),expression(italic("Salvelinus fontinalis"))))+
  xlab("Geographic distance (km)")+
  ylab("Fst/(1-Fst)")+
  theme_bw()+
  facet_wrap(~species,)+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = c(0.2,0.9),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))


ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/figures/comp_isobydist.png",height = 7, width = 12)



ggplot(fst,aes(x = values, y = fst/(1-fst),color = species,fill=species))+
  geom_point(size=3,alpha = 0.8)+
  geom_smooth(method="lm",alpha = 0.3,size = 1,linetype = 2,se = 0.95)+
  scale_color_manual(values = c("#8DA5BEFF","#E49A36FF"),labels = c(expression(italic("Corenus clupeaformis")),expression(italic("Salvelinus fontinalis"))))+
  scale_fill_manual(values = c("#8DA5BEFF","#E49A36FF"),labels = c(expression(italic("Corenus clupeaformis")),expression(italic("Salvelinus fontinalis"))))+
  xlab("Geographic distance (km)")+
  ylab("Fst/(1-Fst)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = c(0.2,0.9),
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/figures/comp_isobydist_merged.png",height = 7, width = 12)


#COCL IBD

#Mantel test

fst_cocl <- fst %>% filter(species == "Coregonus clupeaformis")

#abundance data frame - bray curtis dissimilarity
dist <- as.data.frame.matrix(xtabs(values ~ pop1 + pop2, data=fst_cocl))
dist[dist==0] <- NA

#environmental vector - euclidean distance
mfst <- as.data.frame.matrix(xtabs(fst ~ pop1 + pop2, data=fst_cocl))
mfst[mfst==0] <- NA

mantel(xdis = dist, ydis = mfst, method = "spearman", permutations = 9999, na.rm = TRUE) 

#SAFO IBD

#Mantel test

fst_safo <- fst %>% filter(species == "Salvelinus fontinalis")

#abundance data frame - bray curtis dissimilarity
dist <- as.data.frame.matrix(xtabs(values ~ pop1 + pop2, data=fst_safo))
dist[dist==0] <- NA

#environmental vector - euclidean distance
mfst <- as.data.frame.matrix(xtabs(fst ~ pop1 + pop2, data=fst_safo))
mfst[mfst==0] <- NA

mantel(xdis = dist, ydis = mfst, method = "spearman", permutations = 9999, na.rm = TRUE) 

#

slope_comp <- lm(fst ~ values * species, data = fst)

summary(slope_comp)

####################################
##########  FST heatmap  ###########
####################################

#---> Whitefish

fst_cocl2 <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/fst/fst_cocl.tsv")

colnames(fst_cocl2) <- c("fst","pop1","pop2")

fst_cocl2$pop1 <- gsub("CHI","KAA",fst_cocl2$pop1)
fst_cocl2$pop2 <- gsub("CHI","KAA",fst_cocl2$pop2)

fst_cocl2 %<>% filter(!pop1 == "FIS",!pop2 == "FIS")

fst_cocl2$pop1 <- factor(fst_cocl2$pop1, levels=c("ROG","KAA","BIA","LAG","BEA","SAB","REN","CON","EAS","JAC","RUP","NOT"))
fst_cocl2$pop2 <- factor(fst_cocl2$pop2, levels=c("ROG","KAA","BIA","LAG","BEA","SAB","REN","CON","EAS","JAC","RUP","NOT"))


mat_fst <- as.data.frame.matrix(xtabs(fst~ pop1 + pop2, data=fst_cocl2))

mat_fst[is.na(mat_fst)] <- 0

levels <- c("ROG","KAA","BIA","LAG","BEA","SAB","REN","CON","EAS","JAC","RUP","NOT")

mean(fst_cocl2$fst)
sd(fst_cocl2$fst)


ggplot(fst_cocl2, aes(x = pop1, y = ordered(pop2, levels=rev(levels)), fill= fst)) + 
  geom_tile() +
  geom_text(aes(label = round(fst, 3)),size=8)+
  scale_fill_gradient(low = "white", high = "red", name="Fst") +
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 20, hjust = 1),
        axis.text.y = element_text(size = 20),
        panel.grid = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) + 
  coord_fixed()

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/figures/cocl_fst_heatmap.png",height = 13, width = 13)


#---> Trout

fst_safo2 <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/fst/fst_safo.tsv")

colnames(fst_safo2) <- c("fst","pop1","pop2")

fst_safo2$pop1 <- gsub("CHI","KAA",fst_safo2$pop1)
fst_safo2$pop2 <- gsub("CHI","KAA",fst_safo2$pop2)

fst_safo2 %<>% filter(!pop1 %in% c("FIS","MOA","WAS"),!pop2 %in% c("FIS","MOA","WAS"))

fst_safo2$pop1 <- factor(fst_safo2$pop1, levels=c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","BRO"))
fst_safo2$pop2 <- factor(fst_safo2$pop2, levels=c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","BRO"))


mat_fst <- as.data.frame.matrix(xtabs(fst~ pop1 + pop2, data=fst_safo2))

mat_fst[is.na(mat_fst)] <- 0

levels <- c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","BRO")

mean(fst_safo2$fst)
sd(fst_safo2$fst)


ggplot(fst_safo2, aes(x = pop1, y = ordered(pop2, levels=rev(levels)), fill= fst)) + 
  geom_tile() +
  geom_text(aes(label = round(fst, 3)),size=8)+
  scale_fill_gradient(low = "white", high = "red", name="Fst") +
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 20, hjust = 1),
        axis.text.y = element_text(size = 20),
        panel.grid = element_blank(),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) + 
  coord_fixed()

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/figures/safo_fst_heatmap.png",height = 13, width = 13)




