rm(list = ls())

#####################
#######Library#######
#####################

library(tidyverse)
library(magrittr)
library(paletteer)
library(plotly)
library(vegan)
library(ggridges)

#####################
#########DATA########
#####################

#Import covariance matrix from pcangsd

cov <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/all_maf0.05pctind0.8_maxdepth10_ALL_CHR.singletons.pruned.cov")

#How many PC axis to retain?

pca <- rda(cov)

summary(pca)

#Broken-stick method to find how many PC axis explain more variation than expected by chance

png("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/figure/broken_stick_pca.jpeg")

screeplot(pca)

dev.off()

#Based on the broken stick method, 4 PC axis should be retained.

#Make PCA using covariance matrix

pca <- eigen(cov)

pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

#add column names
nPC<-dim(pca$vectors)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca.mat)<-c(col_PC)

#add rownames
samples <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/ind_list_analysis/ind.list")
rownames(pca.mat)<-samples$V1

#calculate varsum(eigen_mats$values[eigen_mats$values>=0]
var1<-round(pca$values[1]*100/sum(pca$values[pca$values>=0]),2)
var2<-round(pca$values[2]*100/sum(pca$values[pca$values>=0]),2)
var3<-round(pca$values[3]*100/sum(pca$values[pca$values>=0]),2)
var4<-round(pca$values[4]*100/sum(pca$values[pca$values>=0]),2)


#Make dataset to plot
pca_cocl <- data.frame(pca.mat[,1:20])

#sample id column
pca_cocl$id <- samples$V1

#pop column
pca_cocl %<>% mutate(pop = case_when(str_detect(string = id, "RUP") ~ "RUP",
                                 str_detect(string = id, "NOT") ~ "NOT",
                                 str_detect(string = id, "BEA") ~ "BEA",
                                 str_detect(string = id, "FIS") ~ "EAS",
                                 str_detect(string =id, "JAC") ~ "JAC",
                                 str_detect(string = id, "SAB") ~ "SAB",
                                 str_detect(string = id, "BIA") ~ "BIA",
                                 str_detect(string = id, "LAG") & PC2 < 0.05 ~ "LAG",
                                 str_detect(string = id, "LAG") & PC2 > 0.05 ~ "ROG",
                                 str_detect(string = id, "TIL") ~ "REN",
                                 str_detect(string = id, "MOA") ~ "MOA",
                                str_detect(string = id, "ROG") ~ "ROG",
                                str_detect(string = id, "EAS") ~ "EAS",
                                str_detect(string = id, "CHIm") ~ "KAA",
                                str_detect(string = id, "CON") ~ "CON",
                                str_detect(string = id, "BRO") ~ "BRO",
                                str_detect(string = id, "WAS") ~ "WAS",
                                str_detect(string = id, "REN") ~ "REN",
                                str_detect(string = id, "AQU") ~ "AQU"))

pca_cocl %>% group_by(pop) %>% count()




##add info on pop
pop_file <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/pop_sampled.txt",header=T)

pop_file %>% arrange(desc(lat)) 

pca_cocl %<>% left_join(pop_file,by = c("pop" = "pop_code"))

pca_cocl %>% group_by(pop) %>% count()

ggplot(pca_cocl, aes(x = -PC1,y = PC2,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC1 (2.96%)")+
  ylab("PC2 (0.94%)")+
  scale_color_manual(values = rev(c("#A50026FF",
                                "#DD3D2DFF",
                                "#F67E4BFF",
                                "#FDB366FF",
                                "#FEDA8BFF",
                                "#fddbc7",
                                "#EAECCCFF",
                                "#abdda4",
                                "#74add1",
                                "#4575b4",
                                "#313695",
                                "#5e4fa2")))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("figure/pc1_pc2_230k_snp_ld_pruned.jpeg",height=6,width = 8)

ggplot(pca_cocl, aes(x = PC3,y = PC2,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC3 (0.81%)")+
  ylab("PC4 (0.73%)")+
  scale_color_manual(values = rev(c("#A50026FF",
                                    "#DD3D2DFF",
                                    "#F67E4BFF",
                                    "#FDB366FF",
                                    "#FEDA8BFF",
                                    "#fddbc7",
                                    "#EAECCCFF",
                                    "#abdda4",
                                    "#74add1",
                                    "#4575b4",
                                    "#313695",
                                    "#5e4fa2")))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("figure/pc3_pc4_230k_snp_ld_pruned.jpeg",height=6,width = 8)

#######################################
### Check if pattern of missingness ###
#######################################

missing <- read.table("data/missing_data/missing_per_ind.txt")

pca_cocl$missing <- missing$V2

ggplot(pca_cocl, aes(x = -PC1,y = PC2,color = missing))+
  geom_point(alpha=0.4,size = 3)+
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

ggplot(pca_cocl, aes(x = PC1,y = PC3,color = missing))+
  geom_point(alpha=0.8,size = 3)+
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


ggplot(pca_cocl, aes(x = fct_reorder(pop,desc(lat)), y = missing))+
  geom_boxplot()+
  xlab("Population (ordered by latitude)")+
  ylab("Proportion of missing data")+
  theme_bw()

ggsave("figure/missing_data_per_pop.jpeg")

mean(pca_cocl$missing)
mean(pca_cocl$missing) + 2*sd(pca_cocl$missing)

pca_cocl[which(pca_cocl$missing > mean(pca_cocl$missing) + 2*sd(pca_cocl$missing)),]

####################################
#Find best cluster with APCLUSTER
library(apcluster)

df_apclus <- pca_cocl %>% filter(!str_detect(id,"coclNOTs_7023-21|coclRUPs_016-20")) %>% dplyr::select(PC1,PC2,id)

d.apclus <- apcluster(negDistMat(r=1.59),df_apclus,q=0.5)

length(d.apclus@clusters)

cl1a <- cutree(d.apclus, k=8)

plot(cl1a,df_apclus)

# 8 is the best cluster

aggres1a <- aggExCluster(x=d.apclus)

aggres1a@labels <- c("BEA","LAG","KAA","EAS","RUP","NOT","SAB","ROG")

png("figure/cluster_dandogram.png")

plot(aggres1a,showSamples=F)

dev.off()

png("figure/apclus_heatmap_whitefish_8_clus.png")

heatmap(aggres1a,df_apclus)

dev.off()

############################
#Make plot pretty with ggplot
library(data.table)

tidy_apclust_res <- reshape2::melt(d.apclus@clusters)

tidy_apclust_res$value <- as.factor(tidy_apclust_res$value)

df_apclus$rowid <- as.factor(row.names(df_apclus))

plot_apclus <- left_join(df_apclus,tidy_apclust_res,by=c("rowid"="value"))

plot_apclus$L1 <- as.factor(plot_apclus$L1)

plot_apclus %<>% mutate(pop_clus = case_when(L1 == "1" ~ "BEA",
                                           L1 == "2" ~ "LAG",
                                           L1 == "3" ~ "KAA",
                                           L1 == "4" ~ "EAS",
                                           L1 == "5" ~ "RUP",
                                           L1 == "6" ~ "NOT",
                                           L1 == "7" ~ "SAB",
                                           L1 == "8" ~ "ROG"))

sites <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/environmental_data/pop_coord/whitefish_pop_coor_rivers.txt",header=T)

plot_apclus %<>% left_join(sites,by = c("pop_clus" = "SITE_CODE"))

library(plyr)

find_hull <- function(df)df[chull(df$PC1, df$PC2), ]
hulls <- ddply(plot_apclus, "pop_clus", find_hull)

ggplot()+
  geom_point(data = plot_apclus, aes(x=-PC1,y=PC2,colour = reorder(pop_clus,-LAT),fill = reorder(pop_clus,-LAT)),alpha=0.8,size = 3)+
  geom_polygon(data = hulls,aes(x=-PC1,y=PC2,colour = reorder(pop_clus,-LAT),fill = reorder(pop_clus,-LAT)),alpha = 0.5)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  xlab("PC1 (2.96%)")+
  ylab("PC2 (0.94%)")+
  labs(fill='Cluster',color = "Cluster")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_text(size = 18))

ggsave("figure/ap_clus_final.jpeg")

ggplot()+
 geom_polygon(data = hulls,aes(x=-PC1,y=PC2,fill = reorder(pop_clus,-LAT)),alpha = 0.5)+
  geom_point(data = pca_cocl, aes(x=-PC1,y=PC2,colour = fct_reorder(pop,desc(lat))),alpha=0.8,size = 3)+
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  xlab("PC1 (2.96%)")+
  ylab("PC2 (0.94%)")+
  labs(fill='Cluster',color = "Sampling site")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_text(size = 18))

plot_apclus %>% left_join(pca_cocl %>% dplyr::select(id,pop), by="id")


plot_apclus[119,]

pca_cocl %>% filter(id == "coclBIAs_2407-21")

#-------------------------------------------#
#Create a list of sample_id for each cluster

df_apclus %>% filter(PC2 > 0.05 , L1 == "2")


for(i in 1:length(levels(df_apclus$L1))){
  k <- df_apclus %>% filter(L1 == as.character(i)) %>% dplyr::select(id)
  pop_name <- df_apclus %>% filter(L1 == as.character(i)) %>% pull(pop_clus)
  assign(paste0(pop_name[1],".cluster"),k)
  write.table(k, 
              file = paste0("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/clustering_analysis/",pop_name[1],".cluster"),
              quote = F, col.names = F, row.names=F)
}


nrow(KAP.cluster) + nrow(ROG.cluster) + nrow(LAG.cluster) + nrow(BEA.cluster) + nrow(SAB.cluster) + nrow(EAS.cluster) + nrow(RUP.cluster) + nrow(NOT.cluster)


# source pop ref for mixed_stock

write.table(df_apclus %>% dplyr::select(id,pop_clus),"~/Desktop/Baie_James_Paper/02_mixed_stock/whitefish/03_mixed_stock_analysis/source_ref_unit.txt",quote=F,row.names = F)


###

df_apclus %<>% mutate(sampling_site = case_when(str_detect(string = id, "RUP") ~ "RUP",
                                                     str_detect(string = id, "NOT") ~ "NOT",
                                                     str_detect(string = id, "BEA") ~ "BEA",
                                                     str_detect(string = id, "FIS") ~ "EAS",
                                                     str_detect(string =id, "JAC") ~ "JAC",
                                                     str_detect(string = id, "SAB") ~ "SAB",
                                                     str_detect(string = id, "BIA") ~ "BIA",
                                                     str_detect(string = id, "LAG") ~ "LAG",
                                                     str_detect(string = id, "TIL") ~ "REN",
                                                     str_detect(string = id, "MOA") ~ "MOA",
                                                     str_detect(string = id, "ROG") ~ "ROG",
                                                     str_detect(string = id, "EAS") ~ "EAS",
                                                     str_detect(string = id, "CHIm") ~ "KAA",
                                                     str_detect(string = id, "CON") ~ "CON",
                                                     str_detect(string = id, "BRO") ~ "BRO",
                                                     str_detect(string = id, "WAS") ~ "WAS",
                                                     str_detect(string = id, "REN") ~ "REN",
                                                     str_detect(string = id, "AQU") ~ "AQU"))

#-------------------------------#
#Create a list of individuals that assign to one cluster for each sampling site

df_apclus %>% dplyr::count(sampling_site,pop_clus)

#Create a column that specify whether to keep an ind or not based
#on the fact that it assign to a cluster that breads in the specific
#site. We therefore remove most likely F1 migrant from the samples.

site_all_freq <- df_apclus %>% mutate(keep = case_when(sampling_site == "BEA" & pop_clus == "BEA" ~ "keep",
                                                       sampling_site == "BEA" & pop_clus == "BEA" ~ "keep"
                               sampling_site == "CON" & pop_clus == "EAS" ~ "keep",
                               sampling_site == "EAS" & pop_clus == "EAS" ~ "keep",
                               sampling_site == "JAC" & pop_clus == "EAS" ~ "keep",
                               sampling_site == "KAA" & pop_clus == "KAA" ~ "keep",
                               sampling_site == "LAG" & pop_clus == "LAG" ~ "keep",
                               sampling_site == "NOT" & pop_clus == "NOT" ~ "keep",
                               sampling_site == "REN" & pop_clus == "EAS" ~ "keep",
                               sampling_site == "ROG" & pop_clus == "ROG" ~ "keep",
                               sampling_site == "RUP" & pop_clus == "RUP" ~ "keep",
                               sampling_site == "SAB" & pop_clus == "SAB" ~ "keep",
                                T ~ "no_keep")) %>% filter(keep == 'keep')


site_all_freq$sampling_site <- as.factor(site_all_freq$sampling_site)

for(i in 1:length(levels(site_all_freq$sampling_site))){
  k <- site_all_freq %>% filter(sampling_site == as.character(levels(site_all_freq$sampling_site)[i])) %>% dplyr::select(id)
  assign(paste0(levels(site_all_freq$sampling_site)[i],"site.cluster"),k)
  write.table(k, 
              file = paste0("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/clustering_analysis/",levels(site_all_freq$sampling_site)[i],"site.cluster"),
              quote = F, col.names = F, row.names=F)
}



write_rds(plot_apclus,"~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/ap_cluster/ap_cluster_whitefish.Rds")





