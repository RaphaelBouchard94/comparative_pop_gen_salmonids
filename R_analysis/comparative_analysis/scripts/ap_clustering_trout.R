rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(mapmixture)
library(gridExtra)
library(viridis)
library(patchwork)
library(plyr)


#---------------#
#----- Data ----#
#---------------#

#---------------#
#---> Trout <---#
#---------------#

sites <- read.table("data/pop_sampled.txt",header = T)

sites %<>% dplyr::select(pop_code,lat,long) %>% dplyr::rename(site = pop_code)

ind_trout <- read.table("data/ngs_admix/trout/ind.txt")

trout_apclus <- read_rds("data/ap_cluster/ap_cluster_trout.Rds")

trout_apclus %<>% mutate(cluster = case_when(L1 == 1 ~ "SAB",
                                            L1 == 2 ~ "BEA",
                                            L1 == 3 ~ "LAG",
                                            L1 == 4 ~ "AQU",
                                            L1 == 5 ~ "BIA",
                                            L1 == 6 ~ "REN",
                                            L1 == 7 ~ "ROG/KAA",
                                            L1 == 8 ~ "BRO",
                                            L1 == 9 ~ "EAS",
                                            L1 == 10 ~ "CON"))

#---> Make Structure plot with individual clustered based on APclus analysis
#-> K=10

k10_trout <- read.table("data/ngs_admix/trout/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_10_10.qopt")
k10_trout$id <- ind_trout$V1

k10_trout %<>% left_join(trout_apclus %>% dplyr::select(id,cluster), by = "id")

k10_trout %<>% dplyr::select(cluster,id,V1:V10) 


k10_trout_structure <- k10_trout %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V10) %>% 
  group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))


#---> Structure plot K=10

k10_trout_structure$cluster <- ordered(k10_trout_structure$cluster,
                                       levels = c("BRO",
                                                  "EAS",
                                                  "CON",
                                                  "REN",
                                                  "SAB",
                                                  "BEA",
                                                  "AQU",
                                                  "LAG",
                                                  "BIA",
                                                  "ROG/KAA"))




struct_plot_k10_trout <- ggplot(k10_trout_structure, aes(id, prob, fill = pop)) +
  geom_col(width = 2) +
  scale_fill_manual(values = c("#EAECCCFF", #Sabacunica
                               "#FDB366FF", #Conn
                               "#FEDA8BFF", #Renoyer
                               "#DD3D2DFF", #Broadback
                               "#F67E4BFF", #Eastmain
                               "#4A7BB7FF", #Biagawshi
                               "#6EA6CDFF", #La Grande
                               "#98CAE1FF", #Aquatuc
                               "#364B9AFF", #Roggan/Kaapsaoui
                               "#C2E4EFFF")) + #BEAVER
  facet_grid(~cluster,scales = "free_x")+
  ylab("Ancestry")+
  xlab("")+
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        title = element_text(size = 18),
        strip.text = element_text(size = 16),
        legend.position = "none")


struct_plot_k10_trout

#---> APcluster

find_hull <- function(df)df[chull(df$PC1, df$PC2), ]
hulls <- ddply(trout_apclus, "L1", find_hull)

trout_apclus$cluster <- as.factor(trout_apclus$cluster)

#---> Merge APcluster result and admixture 

k10_trout$best_k <- colnames(k10_trout)[apply(k10_trout,1,which.max)]

k10_trout %<>% 
  rowwise() %>% dplyr::mutate(best_k_admix = max(c_across(starts_with("V")), na.rm = TRUE))

k10_trout_to_merge <- k10_trout %>% dplyr::select(id,best_k,best_k_admix)

trout_apclus %<>% left_join(k10_trout_to_merge, by = "id")



trout_apclus %<>% dplyr::mutate(color_per_pop = case_when(cluster == "BRO" ~ "#DD3D2DFF",
                                                         cluster == "EAS" ~ "#F67E4BFF",
                                                         cluster == "CON" ~ "#FDB366FF",
                                                         cluster == "REN" ~ "#FEDA8BFF",
                                                         cluster == "SAB" ~ "#EAECCCFF",
                                                         cluster == "BEA" ~ "#C2E4EFFF",
                                                         cluster == "AQU" ~ "#98CAE1FF",
                                                         cluster == "LAG" ~ "#6EA6CDFF",
                                                         cluster == "BIA" ~ "#4A7BB7FF",
                                                         cluster == "ROG/KAA" ~ "#364B9AFF"))


trout_apclus_plot <- ggplot(data =trout_apclus,aes(x=-PC1,y=PC2, color = cluster, fill = cluster, alpha=best_k_admix))+
  geom_point(size = 3)+
  xlab("PC1")+
  geom_polygon(data = hulls,alpha = 0.5)+
  scale_color_manual(values = c("BRO" = "#DD3D2DFF",
                                "EAS" = "#F67E4BFF",
                                "CON" = "#FDB366FF",
                                "REN" = "#FEDA8BFF",
                                "SAB" = "#EAECCCFF",
                                "BEA" = "#C2E4EFFF",
                                "AQU" = "#98CAE1FF",
                                "LAG" = "#6EA6CDFF",
                                "BIA" = "#4A7BB7FF",
                                "ROG/KAA" = "#364B9AFF"))+
  scale_fill_manual(values = c("BRO" = "#DD3D2DFF",
                               "EAS" = "#F67E4BFF",
                               "CON" = "#FDB366FF",
                               "REN" = "#FEDA8BFF",
                               "SAB" = "#EAECCCFF",
                               "BEA" = "#C2E4EFFF",
                               "AQU" = "#98CAE1FF",
                               "LAG" = "#6EA6CDFF",
                               "BIA" = "#4A7BB7FF",
                               "ROG/KAA" = "#364B9AFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_blank(),
        legend.position = "none")

trout_apclus_plot2 <- ggplot()+
  geom_point(data =trout_apclus,aes(x=-PC1,y=PC2,color = best_k, alpha=best_k_admix), size = 3)+
  scale_fill_manual(values = c("#EAECCCFF", #Sabacunica
                               "#FDB366FF", #Conn
                               "#FEDA8BFF", #Renoyer
                               "#DD3D2DFF", #Broadback
                               "#F67E4BFF", #Eastmain
                               "#4A7BB7FF", #Biagawshi
                               "#6EA6CDFF", #La Grande
                               "#98CAE1FF", #Aquatuc
                               "#364B9AFF", #Roggan/Kaapsaoui
                               "#C2E4EFFF"))+
  scale_color_manual(values = c("#EAECCCFF", #Sabacunica
                                "#FDB366FF", #Conn
                                "#FEDA8BFF", #Renoyer
                                "#DD3D2DFF", #Broadback
                                "#F67E4BFF", #Eastmain
                                "#4A7BB7FF", #Biagawshi
                                "#6EA6CDFF", #La Grande
                                "#98CAE1FF", #Aquatuc
                                "#364B9AFF", #Roggan/Kaapsaoui
                                "#C2E4EFFF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_blank(),
        legend.position = "none")


(trout_apclus_plot+trout_apclus_plot2) / struct_plot_k10_trout

ggsave("figures/apclus_admix_trout.png",width = 12.5)
  
trout_apclus_plot+struct_plot_k10_trout+
  plot_layout(widths = c(.75, 2))

ggsave("figures/apclus_admix_trout2.png",width = 18,height=5)



#Write table
write.table(trout_apclus %>% dplyr::select(id, pop),
            file = "data/wgs_assign/trout/wgs_assign_pop_id.txt",
            quote=F,
            col.names = F,
            row.names = F,
            sep = "\t")








