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


####################################
######## Bathymetric map  ##########
####################################

JB <- getNOAA.bathy(lon1 = -82.084, lon2 = -78,
                        lat1 = 51.021, lat2 = 55, resolution = 1)

#blues <- colorRampPalette(c("red","purple","blue", "cadetblue1","white"))

#plot(JB, image = TRUE,land=TRUE,bpal = list(c(0, max(JB), "grey"),
                                 # c(min(JB),0,blues(100))))
#scaleBathy(JB, deg = 2, x = "bottomleft", inset = 5)
#plot(JB, deep = 0, shallow = 0, step = 0, lwd = 0.4, add = TRUE)
#points(marsites,pch=21,col="yellow",bg=col2alpha("yellow",.9),cex=1.2)

sites <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/fst/pop_sampled_whitefish.txt",header=T)

marsites<- sites %>% select(long,lat) %>% rename(x=long,y=lat)

trans1 <- trans.mat(JB,min.depth = 0)

dist1 <- lc.dist(trans1,marsites,res="dist")

dist1

#out1 <- lc.dist(trans1,marsites,res="path")

#path.profile(out1[[17]],JB,pl=TRUE,
             #main="Path between locations 1 & 3\nProfile with no depth constraint")


#plot(JB, image = TRUE,land=TRUE,bpal = list(c(0, max(JB), "grey"),
                                            c(min(JB),0,blues(100))))
#scaleBathy(JB, deg = 2, x = "bottomleft", inset = 5)
#plot(JB, deep = 0, shallow = 0, step = 0, lwd = 0.4, add = TRUE)
#points(marsites,pch=21,col="red",bg=col2alpha("red",.9),cex=1.2)
#lines(out1[[2]], cex = 0.5, col = "darkorange1")


####################################
##########  FST heatmap  ###########
####################################

fst <- read.csv("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/fst/fst.csv")

fst$pop1 <- factor(fst$pop1, levels=c("KAA","ROG","LAG","BEA","SAB","EAS","RUP","NOT"))
fst$pop2 <- factor(fst$pop2, levels=c("KAA","ROG","LAG","BEA","SAB","EAS","RUP","NOT"))


mat_fst <- as.data.frame.matrix(xtabs(fst~ pop1 + pop2, data=fst))

mat_fst[is.na(mat_fst)] <- 0

levels <-c("KAA","ROG","LAG","BEA","SAB","EAS","RUP","NOT")

mean(fst$fst)
sd(fst$fst)


ggplot(fst, aes(x = pop1, y = ordered(pop2, levels=rev(levels)), fill= fst)) + 
  geom_tile() +
  geom_text(aes(label = round(fst, 3)),size=8)+
  annotate("text",x=5.5,y=7,label="mean Fst = 0.055±0.030",size=8)+
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

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/figure/fst_heatmap.jpeg",height = 10, width = 10)


####################################
###  FST Isolation by distance  ####
####################################

mat_dist<-as.matrix(dist1)

dfdist <- setNames(melt(mat_dist), c('pop1', 'pop2', 'values'))

str(dfdist)

dfdist %<>% mutate(pop1_code = case_when(pop1 == 1 ~ "KAA",
                                        pop1 == 2 ~ "ROG",
                                        pop1 == 3 ~ "LAG",
                                        pop1 == 4 ~ "BEA",
                                        pop1 == 5 ~ "SAB",
                                        pop1 == 6 ~ "EAS",
                                        pop1 == 7 ~ "RUP",
                                        pop1 == 8 ~ "NOT"),
                  pop2_code = case_when(pop2 == 1 ~ "KAA",
                                        pop2 == 2 ~ "ROG",
                                        pop2 == 3 ~ "LAG",
                                        pop2 == 4 ~ "BEA",
                                        pop2 == 5 ~ "SAB",
                                        pop2 == 6 ~ "EAS",
                                        pop2 == 7 ~ "RUP",
                                        pop2 == 8 ~ "NOT")) %>% select(-(pop1:pop2))

fst %<>% left_join(dfdist,by=c("pop1"="pop1_code","pop2"="pop2_code"))

#Mantel test

#abundance data frame - bray curtis dissimilarity
dist <- as.data.frame.matrix(xtabs(values ~ pop1 + pop2, data=fst))
dist[dist==0] <- NA

#environmental vector - euclidean distance
mfst <- as.data.frame.matrix(xtabs(fst ~ pop1 + pop2, data=fst))
mfst[mfst==0] <- NA

mantel(xdis = dist, ydis = mfst, method = "spearman", permutations = 9999, na.rm = TRUE) 


#Isolation by distance graph

fst %<>% mutate(label_comp = paste(pop1,pop2,sep="_"))

ggplot(fst,aes(x = values, y = fst/(1-fst)))+
  geom_point()+
  geom_smooth(method="lm",color="black")+
  annotate("text",x = 100, y = 0.14,label="Mantel r : 0.86",size=6)+
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


ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/figure/isobydist.jpeg",height = 7, width = 7)







