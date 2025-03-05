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

####################################
######## Bathymetric map  ##########
####################################

JB <- getNOAA.bathy(lon1 = -82.084, lon2 = -78,
                    lat1 = 51.021, lat2 = 55, resolution = 1)

#blues <- colorRampPalette(c("red","purple","blue", "cadetblue1","white"))

#plot(JB, image = TRUE,land=TRUE,bpal = list(c(0, max(JB), "grey"),
                                            c(min(JB),0,blues(100))))
#scaleBathy(JB, deg = 2, x = "bottomleft", inset = 5)
#plot(JB, deep = 0, shallow = 0, step = 0, lwd = 0.4, add = TRUE)
#points(marsites,pch=21,col="yellow",bg=col2alpha("yellow",.9),cex=1.2)

sites <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/fst/pop_sampled.txt",header=T)

marsites<- sites %>% select(long,lat) %>% rename(x=long,y=lat)

trans1 <- trans.mat(JB,min.depth = 0)

dist1 <- lc.dist(trans1,marsites,res="dist")
dist1

#out1 <- lc.dist(trans1,marsites,res="path")

#path.profile(out1[[28]],JB,pl=TRUE,
             #main="Path between locations 1 & 3\nProfile with no depth constraint")

####################################
##########  FST heatmap  ###########
####################################

fst <- read.csv("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/fst/fst_trout.csv")

fst$pop1 <- factor(fst$pop1, levels=c("KAA","BIA","LAG","AQU","BEA","SAB","REN","EAS","BRO"))
fst$pop2 <- factor(fst$pop2, levels=c("KAA","BIA","LAG","AQU","BEA","SAB","REN","EAS","BRO"))


mat_fst <- as.data.frame.matrix(xtabs(fst~ pop1 + pop2, data=fst))

mat_fst[is.na(mat_fst)] <- 0

levels <-c("KAA","BIA","LAG","AQU","BEA","SAB","REN","EAS","BRO")

mean(fst$fst)
sd(fst$fst)


ggplot(fst, aes(x = pop1, y = ordered(pop2, levels=rev(levels)), fill= fst)) + 
  geom_tile() +
  geom_text(aes(label = round(fst, 3)),size = 8)+
  annotate("text",x=5.5,y=7,label="mean Fst = 0.077±0.039",size=8)+
  scale_fill_gradient(low = "white", high = "red", name="Fst") +
  theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 18, hjust = 1),
        axis.text.y = element_text(size = 18),
        panel.grid = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18)) + 
  coord_fixed()

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/fst_heatmap.jpeg",height = 10, width = 10)


####################################
###  FST Isolation by distance  ####
####################################

mat_dist<-as.matrix(dist1)

dfdist <- setNames(melt(mat_dist), c('pop1', 'pop2', 'values'))

dfdist %<>% mutate(pop1_code = case_when(pop1 == 1 ~ "KAA",
                                         pop1 == 2 ~ "BIA",
                                         pop1 == 3 ~ "LAG",
                                         pop1 == 4 ~ "AQU",
                                         pop1 == 5 ~ "BEA",
                                         pop1 == 6 ~ "SAB",
                                         pop1 == 7 ~ "REN",
                                         pop1 == 8 ~ "EAS",
                                         pop1 == 9 ~ "BRO"),
                   pop2_code = case_when(pop2 == 1 ~ "KAA",
                                         pop2 == 2 ~ "BIA",
                                         pop2 == 3 ~ "LAG",
                                         pop2 == 4 ~ "AQU",
                                         pop2 == 5 ~ "BEA",
                                         pop2 == 6 ~ "SAB",
                                         pop2 == 7 ~ "REN",
                                         pop2 == 8 ~ "EAS",
                                         pop2 == 9 ~ "BRO")) %>% select(-(pop1:pop2))

fst %<>% left_join(dfdist,by=c("pop1"="pop1_code","pop2"="pop2_code"))

fst %<>% mutate(comp_tag = paste(pop1,pop2,sep="_"))

#Mantel test

#abundance data frame - bray curtis dissimilarity
dist <- as.data.frame.matrix(xtabs(values ~ pop1 + pop2, data=fst))
dist[dist==0] <- NA

#environmental vector - euclidean distance
mfst <- as.data.frame.matrix(xtabs(fst ~ pop1 + pop2, data=fst))
mfst[mfst==0] <- NA

mantel(xdis = dist, ydis = mfst, method = "spearman", permutations = 9999, na.rm = TRUE) 


fst %>% 
  ggplot(aes(x = values, y = fst/(1-fst)))+
  geom_point()+
  geom_text(aes(label=comp_tag))+
  geom_smooth(method="lm")+
  xlab("Distance (km)")+
  ylab("Fst/(1-Fst)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/isobydist.jpeg",height = 7, width = 7)




