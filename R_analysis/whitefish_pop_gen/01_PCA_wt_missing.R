rm(list = ls())

#####################
#######Library#######
#####################

library(tidyverse)
library(magrittr)
library(paletteer)
library(plotly)
library(vegan)


#####################
#########DATA########
#####################

#Import covariance matrix from pcangsd

cov <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/pca_without_missing/all_maf0.05pctind0.8_maxdepth10_ALL_CHR.singletons.pruned.cov")

#How many PC axis to retain?

pca <- rda(cov)

screeplot(pca)

summary(pca)

#Broken-stick method to find how many PC axis explain more variation than expected by chance

#png("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/figure/broken_stick_pca.jpeg",res=70)

screeplot(pca,bstick = T)

#dev.off()

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
samples <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/pca_without_missing/bam2keep.filelist")
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
pop_file <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/pop_sampled.txt",header=T)

pop_file %>% arrange(desc(lat)) 

pca_cocl %<>% left_join(pop_file,by = c("pop" = "pop_code"))

pca_cocl %>% group_by(pop) %>% count()

ggplot(pca_cocl, aes(x = -PC1,y = PC2,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC1 (2.96%)")+
  ylab("PC2 (0.94%)")+
  scale_color_manual(values = c("#364B9AFF",
                                "#4A7BB7FF",
                                "#6EA6CDFF",
                                "#98CAE1FF",
                                "#C2E4EFFF",
                                "#EAECCCFF",
                                "#FDDBC7FF",
                                "#FEDA8BFF",
                                "#FDB366FF",
                                "#F67E4BFF",
                                "#DD3D2DFF",
                                "#A50026FF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("figure/pc1_pc2_230k_snp_ld_pruned.jpeg",height=6,width = 8)


ggplot(pca_cocl, aes(x = PC3,y = PC4,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC3 (0.83%)")+
  ylab("PC4 (0.76%)")+
  scale_color_manual(values = c("#364B9AFF",
                                "#4A7BB7FF",
                                "#6EA6CDFF",
                                "#98CAE1FF",
                                "#C2E4EFFF",
                                "#EAECCCFF",
                                "#FDDBC7FF",
                                "#FEDA8BFF",
                                "#FDB366FF",
                                "#F67E4BFF",
                                "#DD3D2DFF",
                                "#A50026FF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("figure/pc3_pc4_230k_snp_ld_pruned.jpeg",height=6,width = 8)

ggplot(pca_cocl, aes(x = PC5,y = PC6,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC3 (0.81%)")+
  ylab("PC4 (0.73%)")+
  scale_color_viridis_d()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggplot(pca_cocl, aes(x = PC7,y = PC8,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC3 (0.81%)")+
  ylab("PC4 (0.73%)")+
  scale_color_viridis_d()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


#weird ind for clustering
pca_cocl %>% filter(between(PC1,-0.1,0) & between(PC2,0.05,0.12))

pca_cocl %>% 
ggplot(aes(x = PC1,y = PC2,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC1 (2.96%)")+
  ylab("PC2 (0.94%)")+
  scale_color_viridis_d()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

####################################
#Find best cluster with APCLUSTER
library(apcluster)

df_apclus <- pca_cocl %>% filter(!str_detect(id,"coclNOTs_7023-21|coclRUPs_016-20")) %>% dplyr::select(PC1,PC2)

d.apclus <- apcluster(negDistMat(df_apclus,r=1.9))

length(d.apclus@clusters)

cl1a <- cutree(d.apclus, k=8)

plot(cl1a,df_apclus,
     col=viridis::viridis(9),
     bg=viridis::viridis(9))

############################
#Make plot pretty with ggplot
library(data.table)
library(plyr)
  
deme_df$clus <- pam_res$clustering

deme_df$clus <- as.factor(deme_df$clus)

deme_df %<>% mutate(pop_clus = case_when(clus == "1" ~ "BEA",
                                           clus == "2" ~ "SAB",
                                           clus == "3" ~ "LAG",
                                           clus == "4" ~ "KAA",
                                           clus == "5" ~ "ROG",
                                           clus == "6" ~ "RUP",
                                           clus == "7" ~ "EAS",
                                           clus == "8" ~ "NOT"))

find_hull <- function(df) df[chull(df$PC1, df$PC2), ]

hulls <- ddply(deme_df, "pop_clus", find_hull)

hulls %<>% left_join(pop_file, by=c("pop_clus"="pop_code")) 

deme_df %<>% left_join(pop_file, by=c("pop_clus"="pop_code"))

ggplot(deme_df,aes(x=-PC1,y=PC2,fill = fct_reorder(pop_clus,-lat),color = fct_reorder(pop_clus,-lat)))+
  geom_point(alpha=0.8,size = 3)+
  geom_polygon(data = hulls,alpha = 0.5)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_blank())

ggsave("figure/pam_clus.jpeg")

#--------------------#

deme_df$origin <- pca_cocl$pop

write.table(deme_df,"data/individual_assignment_analysis/ind_assignment.txt")

deme_df_prop <- deme_df %>% dplyr::group_by(pop_clus) %>% dplyr::summarise(sum_prop = sum(origin))


ggplot(deme_df,aes(x=fct_reorder(pop_clus,-lat),fill=origin))+
  geom_bar()



