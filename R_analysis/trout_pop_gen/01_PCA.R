rm(list = ls())

#####################
#######Library#######
#####################

library(tidyverse)
library(magrittr)
library(paletteer)
library(plotly)
library(vegan)
library(ade4)
library(MASS)


#####################
#########DATA########
#####################

cov <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/pcangsd/all_maf0.05pctind0.85_maxdepth10_ALL_CHR_singletons.pruned.cov")

#pca <- rda(cov)

#screeplot(pca,bstick = T)

#Make PCA using covariance matrix
pca <- eigen(cov)
pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

#add column names
nPC<-dim(pca$vectors)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca.mat)<-c(col_PC)

#add rownames
samples <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/pcangsd/ind.list")
rownames(pca.mat)<-samples$V1

#calculate varsum(eigen_mats$values[eigen_mats$values>=0]
var1<-round(pca$values[1]*100/sum(pca$values[pca$values>=0]),2)
var2<-round(pca$values[2]*100/sum(pca$values[pca$values>=0]),2)
var3<-round(pca$values[3]*100/sum(pca$values[pca$values>=0]),2)
var4<-round(pca$values[4]*100/sum(pca$values[pca$values>=0]),2)


#Make dataset to plot
pca_trout <- data.frame(pca.mat[,1:20])

#sample id column
pca_trout$id <- samples$V1

#pop column
pca_trout %<>% mutate(pop = case_when(str_detect(string = id, "BEA") ~ "BEA",
                                     str_detect(string = id, "FIS") ~ "EAS",
                                     str_detect(string = id, "JAC") ~ "JAC",
                                     str_detect(string = id, "SAB") ~ "SAB",
                                     str_detect(string = id, "BIA") ~ "BIA",
                                     str_detect(string = id, "LAG") ~ "LAG",
                                     str_detect(string = id, "TIL") ~ "TIL",
                                     str_detect(string = id, "MOA") ~ "SAB",
                                     str_detect(string = id, "ROG") ~ "ROG",
                                     str_detect(string = id, "EAS") ~ "EAS",
                                     str_detect(string = id, "CHIm") ~ "KAA",
                                     str_detect(string = id, "CON") ~ "CON",
                                     str_detect(string = id, "AQU") ~ "AQU",
                                     str_detect(string = id, "REN") ~ "REN",
                                     str_detect(string = id, "WAS") ~ "BRO",
                                     str_detect(string = id, "BRO") ~ "BRO"))

##add info on pop
pop_file <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/pop_sampled.txt",header=T)

pop_file %>% arrange(desc(lat)) 

pca_trout %<>% left_join(pop_file,by = c("pop" = "pop_code"))

pca_trout %>% group_by(pop) %>% count()

ggplot(pca_trout, aes(x = -PC1,y = PC2,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC1 (1.93%)")+
  ylab("PC2 (1.42%)")+
  scale_color_manual("Cluster",values = c("#C62828FF",
                                          "#F44336FF",
                                          "darkorchid2",#BIA
                                          "#673AB7FF",#LAG
                                          "lightskyblue2",#AQU
                                          "#3F51B5FF",#BEA
                                          "#2196F3FF",#SAB
                                          "darkgreen",#REN
                                          "#8BC34AFF",#CON
                                          "#4CAF50FF",#EAS
                                          "lightgoldenrod2"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/pc1_pc2_118k_snp_ld_pruned.jpeg",height=6,width = 8)

ggplot(pca_trout, aes(x = PC3,y = PC4,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC3 (1.21%)")+
  ylab("PC4 (1.05%)")+
  scale_color_manual("Cluster",values = c("#C62828FF",
                                          "#F44336FF",
                                          "darkorchid2",#BIA
                                          "#673AB7FF",#LAG
                                          "lightskyblue2",#AQU
                                          "#3F51B5FF",#BEA
                                          "#2196F3FF",#SAB
                                          "darkgreen",#REN
                                          "#8BC34AFF",#CON
                                          "#4CAF50FF",#EAS
                                          "lightgoldenrod2"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/pc3_pc4_118k_snp_ld_pruned.jpeg",height=6,width = 8)

#######################################
### Check if pattern of missingness ###
#######################################

missing <- read.table("data/missing_data/missing_per_ind.txt")

pca_trout$missing <- missing$V2

ggplot(pca_trout, aes(x = -PC1,y = PC2,color = missing))+
  geom_point(alpha=1,size = 3)+
  xlab("PC1 (2.96%)")+
  ylab("PC2 (0.94%)")+
  scale_color_viridis_c()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggplot(pca_trout, aes(x = PC3,y = PC4,color = missing))+
  geom_point(alpha=1,size = 3)+
  xlab("PC1 (2.96%)")+
  ylab("PC2 (0.94%)")+
  scale_color_viridis_c()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggplot(pca_trout, aes(x = fct_reorder(pop,desc(lat)), y = missing))+
  geom_boxplot()+
  xlab("Population (ordered by latitude)")+
  ylab("Proportion of missing data")+
  theme_bw()

ggsave("figures/missing_data_per_pop.jpeg")

mean(pca_trout$missing)
mean(pca_trout$missing) + 2*sd(pca_trout$missing)

nrow(pca_trout[which(pca_trout$missing > mean(pca_trout$missing) + 2*sd(pca_trout$missing)),])

ggplot(pca_trout,aes(x=missing))+
  geom_histogram(color="black",bins=50)


###########################
#####K-mean clustering#####
###########################
library(clValid)
library(factoextra)
library(dbscan)
library(fpc)

deme_df_ap<-pca_trout %>% dplyr::select(PC1,PC2,PC3,PC4)

#######
#What is the best number of cluster?
library(clusterSim)
library(cluster)

set.seed(7)
gap <- clusGap(deme_df_ap, pam, 14, B = 500,verbose = interactive(),d=1)

fviz_nbclust(df_apclus, pam, method = "silhouette", k.max = 14)

dfgap <- as.data.frame(gap$Tab)

#write.table(dfgap,"~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/clustering_analysis/gap_stats.txt",quote=F,row.names=F)

dfgap <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/clustering_analysis/gap_stats.txt",header=T)

dfgap$k <- as.factor(c(1:14))

ggplot(dfgap,aes(x=k,y=gap))+
  geom_point()+
  geom_errorbar(aes(ymin=gap-SE.sim,ymax=gap+SE.sim),width=0.3)+
  geom_vline(xintercept = which.max(dfgap$gap+dfgap$SE.sim), lty=5)+
  ylab("GAP statistic")+
  xlab("Number of cluster K")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/gapstat_ncluster.png",
       width = 6,height = 5)

####################################
#Find best cluster with APCLUSTER
library(apcluster)

df_apclus <- pca_trout %>% dplyr::select(PC1,PC2,PC3,id)

d.apclus <- apclusterK(negDistMat(), df_apclus[,-4],9)

length(d.apclus@clusters)

# 8 is the best cluster

cl1a <- cutree(d.apclus, k=9)

png("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/apcluster.png",res=70)

plot(cl1a,df_apclus,
     col=c("#3F51B5FF",
           "#FF9800FF",
           "#FFEB3BFF",
           "#C62828FF",
           "#29B6F6FF",
           "#673AB7FF",
           "#4CAF50FF",
           "#F44336FF"),
     bg=c("#3F51B5FF",
          "#FF9800FF",
          "#FFEB3BFF",
          "#C62828FF",
          "#29B6F6FF",
          "#673AB7FF",
          "#4CAF50FF",
          "#F44336FF"))

dev.off()


#png("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/apcluster.png",res=70)

#heatmap(d.apclus)

#dev.off()

############################
#Make plot pretty with ggplot
library(data.table)

tidy_apclust_res <- reshape2::melt(d.apclus@clusters)

tidy_apclust_res$value <- as.factor(tidy_apclust_res$value)

df_apclus$rown <- as.factor(row.names(df_apclus))

df_apclus %<>% left_join(tidy_apclust_res,by=c("rown"="value"))

df_apclus$L1 <- as.factor(df_apclus$L1)


df_apclus %<>% mutate(pop_clus = case_when(L1 == "1" ~ "REN",
                                           L1 == "2" ~ "SAB",
                                           L1 == "3" ~ "AQU",
                                           L1 == "4" ~ "BIA",
                                           L1 == "5" ~ "LAG",
                                           L1 == "6" ~ "BEA",
                                           L1 == "7" ~ "KAA",
                                           L1 == "8" ~ "EAS",
                                           L1 == "9" ~ "BRO"))

sites <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/fst/pop_sampled.txt",header=T)


df_apclus %<>% left_join(sites,by = c("pop_clus" = "pop_code"))

library(plyr)

find_hull <- function(df)df[chull(df$PC1, df$PC2), ]
hulls <- ddply(df_apclus, "pop_clus", find_hull)

ggplot(df_apclus,aes(x=-PC1,y=PC2,colour = reorder(pop_clus,-lat),fill = reorder(pop_clus,-lat)))+
  geom_point(alpha=0.8,size = 3)+
  geom_polygon(data = hulls,alpha = 0.5)+
  xlab("PC1")+
  scale_color_manual("Cluster",values = c("#C62828FF",#KAA
                                "darkorchid2",#BIA
                                "#673AB7FF",#LAG
                                "lightskyblue2",#AQU
                                "#3F51B5FF",#BEA
                                "#2196F3FF",#SAB
                                "darkgreen",#REN
                                "#4CAF50FF",#EAS
                                "lightgoldenrod2"))+#BRO
  scale_fill_manual("Cluster",values = c("#C62828FF",#KAA
                               "darkorchid2",#BIA
                               "#673AB7FF",#LAG
                               "lightskyblue2",#AQU
                               "#3F51B5FF",#BEA
                               "#2196F3FF",#SAB
                               "darkgreen",#REN
                               "#4CAF50FF",#EAS
                               "lightgoldenrod2"))+#BRO
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_text(size = 18))

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/ap_cluster_res_k9.png",height=6,width = 8)


#

for(i in 1:length(levels(df_apclus$L1))){
  k <- df_apclus %>% filter(L1 == as.character(i)) %>% dplyr::select(id)
  pop_name <- df_apclus %>% filter(L1 == as.character(i)) %>% pull(pop_clus)
  assign(paste0(pop_name[1],".cluster"),k)
  write.table(k, 
              file = paste0("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/clustering_analysis/",pop_name[1],".cluster"),
              quote = F, col.names = F, row.names=F)
}


#nrow(AQU.cluster) + nrow(BEA.cluster) + nrow(BIA.cluster) + nrow(BRO.cluster) + nrow(EAS.cluster) + nrow(KAP.cluster) + nrow(LAG.cluster) + nrow(REN.cluster) + nrow(SAB.cluster)
#Total = 374


