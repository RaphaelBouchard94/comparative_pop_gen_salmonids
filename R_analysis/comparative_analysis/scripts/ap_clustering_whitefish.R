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

#-------------------#
#---> Whitefish <---#
#-------------------#

sites <- read.table("data/pop_sampled.txt",header = T)

sites %<>% dplyr::select(pop_code,lat,long) %>% dplyr::rename(site = pop_code)

ind_whitefish <- read.table("data/ngs_admix/whitefish/ind2.txt")

whitefish_apclus <- read_rds("data/ap_cluster/ap_cluster_whitefish.Rds")

whitefish_apclus %<>% mutate(cluster = case_when(L1 == 1 ~ "BEA",
                                             L1 == 2 ~ "LAG",
                                             L1 == 3 ~ "KAA",
                                             L1 == 4 ~ "EAS-C",
                                             L1 == 5 ~ "RUP",
                                             L1 == 6 ~ "NOT",
                                             L1 == 7 ~ "SAB",
                                             L1 == 8 ~ "ROG"))

#---> Make Structure plot with individual clustered based on APclus analysis
#-> K=8

#18

k8_whitefish <- read.table("data/ngs_admix/whitefish/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_8_18.qopt")
k8_whitefish$id <- ind_whitefish$V1

k8_whitefish %<>% left_join(whitefish_apclus %>% dplyr::select(id,cluster), by = "id")

k8_whitefish %<>% dplyr::select(cluster,id,V1:V8) 


k8_whitefish_structure <- k8_whitefish %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V8) %>% 
  group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
                assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))


#---> Structure plot K=8

k8_whitefish_structure$cluster <- ordered(k8_whitefish_structure$cluster,
                                       levels = c("NOT",
                                                  "RUP",
                                                  "EAS-C",
                                                  "SAB",
                                                  "BEA",
                                                  "LAG",
                                                  "KAA",
                                                  "ROG"))

struct_plot_k8_whitefish <- ggplot(k8_whitefish_structure %>% filter(!is.na(cluster)), aes(id, prob, fill = pop)) +
  geom_col(width = 2) +
  scale_fill_manual(values = c("#C2E4EFFF", #BEAVER
                               "#F67E4BFF", #RUPERT
                               "#364B9AFF", #ROGGAN
                               "#A50026FF", #NOTTAWAY
                               "#4575B4FF", #KAAPSAOUI
                               "#EAECCCFF", #SABACUNICA
                               "#6EA6CDFF", #LAGRANDE
                               "#FDB366FF"))+ #EASTMAIN
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

struct_plot_k8_whitefish 

#---> APcluster


#---> Merge APcluster result and admixture 

k8_whitefish %<>%
  rowwise() %>%
  dplyr::mutate(best_k = paste0(names(.)[c_across(starts_with("V")) == max(c_across(starts_with("V")))], collapse = '_'))

k8_whitefish %<>% 
  rowwise() %>% dplyr::mutate(best_k_admix = max(c_across(starts_with("V")), na.rm = TRUE))

k8_whitefish_to_merge <- k8_whitefish%>% dplyr::select(id,best_k,best_k_admix)

k8_whitefish_to_merge %<>% 
  group_by(id) %>% distinct()

whitefish_apclus %<>% left_join(k8_whitefish_to_merge, by = "id")

find_hull <- function(df)df[chull(df$PC1, df$PC2), ]
hulls <- ddply(whitefish_apclus %>% filter(!is.na(best_k)), "L1", find_hull)

whitefish_apclus$cluster <- as.factor(whitefish_apclus$cluster)


pal_custom3 <- c("#C2E4EFFF",#BEAVER
                 "#6EA6CDFF", #LAGRANDE
                 "#4575B4FF", #KAPSAOUI
                 "#FDB366FF", #EASTMAIN
                 "#F67E4BFF", #RUPERT
                 "#A50026FF", #NOTTAWAY
                 "#EAECCCFF", #SABACUNICA
                 "#364B9AFF") #ROGGAN





whitefish_apclus_plot <- ggplot(data =whitefish_apclus %>% filter(!is.na(best_k)),aes(x=-PC1,y=PC2,color = L1, fill = L1, alpha=best_k_admix))+
  geom_point(size = 3)+
  xlab("PC1")+
  geom_polygon(data = hulls,alpha = 0.5)+
  scale_fill_manual(values = pal_custom3)+
  scale_color_manual(values = pal_custom3)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_blank(),
        legend.position = "none")

whitefish_apclus_plot2 <-ggplot()+
  geom_point(data =whitefish_apclus %>% filter(!is.na(best_k)),aes(x=-PC1,y=PC2,color = best_k, alpha=best_k_admix), size = 3)+
  scale_fill_manual(values = c("#C2E4EFFF", #BEAVER
                               "#F67E4BFF", #RUPERT
                               "#364B9AFF", #ROGGAN
                               "#A50026FF", #NOTTAWAY
                               "#4575B4FF", #KAAPSAOUI
                               "#EAECCCFF", #SABACUNICA
                               "#6EA6CDFF", #LAGRANDE
                               "#FDB366FF"))+
  scale_color_manual(values = c("#C2E4EFFF", #BEAVER
                                "#F67E4BFF", #RUPERT
                                "#364B9AFF", #ROGGAN
                                "#A50026FF", #NOTTAWAY
                                "#4575B4FF", #KAAPSAOUI
                                "#EAECCCFF", #SABACUNICA
                                "#6EA6CDFF", #LAGRANDE
                                "#FDB366FF"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_blank(),
        legend.position = "none")


(whitefish_apclus_plot+whitefish_apclus_plot2) / struct_plot_k8_whitefish 

ggsave("figures/apclus_admix_whitefish.png",width = 12.5)


whitefish_apclus_plot + struct_plot_k8_whitefish +
  plot_layout(widths = c(.75, 2))

ggsave("figures/apclus_admix_whitefish2.png",width = 18,height=5)


whitefish_apclus_tmp <- whitefish_apclus %>% 
  dplyr::group_by(id) %>% 
  dplyr::select(id,cluster) %>% 
  distinct() 
    
    
wgs_assign_file <- ind_whitefish %>% left_join(whitefish_apclus_tmp,by=c("V1"="id"))


wgs_assign_file %<>% mutate(cluster2 = case_when(is.na(cluster) & str_detect(V1,"RUP") ~ "Rupert",
                                                is.na(cluster) & str_detect(V1,"NOT") ~ "Nottaway",
                                                T ~ cluster))

write.table(wgs_assign_file %>% dplyr::select(V1,cluster2),
            file = "data/wgs_assign/whitefish/whitefish_wgs_assign_pop_id.txt",
            quote=F,
            col.names = F,
            row.names = F,
            sep = "\t")





