rm(list = ls())

####################################
###########Library##################
####################################

library(tidyverse)
library(magrittr)
library(data.table)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(ggrepel)

####################################
######  PCA  Baie James########
####################################
pca <- read.table("~/Desktop/GT-seq/GT_seq_BT/james_bay/all_maf0.05_pctind0.65_maxdepth10.cov.pca")

jb <- pca

jb$id <- row.names(pca)


jb %<>% mutate(id_new = gsub("/project/lbernatchez/01_projects/Projects/Raphael_Bouchard/2022-05-10_wgs_sample_preparation_safo/09_no_overlap","",jb$id))
jb %<>% mutate(id = gsub("_1.trimmed.sorted.bam.dedup.bam.no_overlap.bam","",jb$id_new))
jb %<>% mutate(id_new2 = gsub("/NS.1834.002.N...........","",jb$id))

jb %<>% dplyr::select(PC1:PC4,id_new2) %>% 
  separate(id_new2, into = c("pop","id"), sep = "_")

jb %<>% 
  mutate(pop_new = str_sub(pop,-4,-2))



#Plot the PCA

jb %>% 
  ggplot(aes(x = PC1,y = PC2,col= pop_new))+
  geom_point(size=3)+
  geom_text(data=subset(jb, PC1 < 0.11 & PC2 < 0.03 & PC2 > -0.05),
            aes(x = PC1,y = PC2,label=id),
            nudge_x = 0.05)+
  xlab("PC1 (8.72%)")+
  ylab("PC2 (1.79%)")+
  scale_colour_manual(values=c("#e31a1c","#1f78b4","forestgreen"))+
  theme_classic()+
  theme(legend.title = element_blank())+
  theme(title = element_text(size = 28),
        axis.text=element_text(size=26),
        axis.title = element_text(size=28),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.background = element_blank(),
        strip.text.x = element_text(size = 26),
        panel.border = element_rect(colour = "black", fill = NA))







jb %>% 
  ggplot(aes(x = PC2,y = PC3,col= pop_new))+
  geom_point(size=3)+
  xlab("PC1 (8.72%)")+
  ylab("PC2 (1.79%)")+
  scale_colour_manual(values=c("#e31a1c","#1f78b4","forestgreen"))+
  theme_classic()+
  theme(legend.title = element_blank())+
  theme(title = element_text(size = 28),
        axis.text=element_text(size=26),
        axis.title = element_text(size=28),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        legend.background = element_blank(),
        strip.text.x = element_text(size = 26),
        panel.border = element_rect(colour = "black", fill = NA))






b

a+b + plot_annotation(tag_levels = 'A')

ggsave("~/Desktop/james_bay_pca_pc1_pc2.jpeg",width = 16,height = 8)

