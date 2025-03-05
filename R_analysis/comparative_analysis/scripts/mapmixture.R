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
library(sf)
library(ggsn)

#---------------#
#----- Data ----#
#---------------#

#---> Format sites data
sites <- read.table("data/pop_sampled.txt",header = T)

sites %<>% dplyr::select(pop_code,lat,long) %>% dplyr::rename(site = pop_code)

ind_trout <- read.table("data/ngs_admix/trout/ind.txt")
ind_whitefish <- read.table("data/ngs_admix/whitefish/ind.txt")

#---> Load river data

rog <- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/roggan/")
kaa <- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/kapsaoui/")
bia <- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/bia/")
lagr <- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/lagrande/")
aqu<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/aquatuc/")
bea<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/beaver/")
sab<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/sabacunica/")
sabl<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/sabacunica_l/")
ren<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/renoyer/")
conn<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/conn/")
eas<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/eastmain/")
jac<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/jack/")
rup<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/rupert/")
bro<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/broadback/")
not<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/nottaway/")

#---------------#
#---> Trout <---#
#---------------#

#---> Start by making dataframes for the relevant structure plot 

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

trout_apclus %<>% dplyr::select(id,cluster)

#-> K=2

k2_trout <- read.table("data/ngs_admix/trout/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_2_18.qopt")
k2_trout$ind <- ind_trout$V1

k2_trout %<>% left_join(trout_apclus,by=c("ind"="id"))

k2_trout %<>% mutate(site = case_when(str_detect(ind,"ROG") ~ "ROG",
                                     str_detect(ind,"CHIm") ~ "KAA",
                                     str_detect(ind,"BIA") ~ "BIA",
                                     str_detect(ind,"LAG") ~ "LAG",
                                     str_detect(ind,"AQU") ~ "AQU",
                                     str_detect(ind,"BEA") ~ "BEA",
                                     str_detect(ind,"SAB") ~ "SAB",
                                     str_detect(ind,"MOA") ~ "SAB",
                                     str_detect(ind,"REN") ~ "REN",
                                     str_detect(ind,"CON") ~ "CON",
                                     str_detect(ind,"EAS") ~ "EAS",
                                     str_detect(ind,"FIS") ~ "EAS",
                                     str_detect(ind,"WASm") & cluster == "BRO" ~ "BRO",
                                     str_detect(ind,"WASm") & cluster == "EAS" ~ "EAS",
                                     str_detect(ind,"WASm") & cluster == "CON" ~ "CON",
                                     str_detect(ind,"WASm") & cluster == "REN" ~ "REN",
                                     str_detect(ind,"RUP") ~ "RUP",
                                     str_detect(ind,"BRO") ~ "BRO",
                                     str_detect(ind,"NOT") ~ "NOT"))


k2_trout %>% group_by(site) %>% dplyr::count()

k2_trout %<>% dplyr::select(site,ind,V1,V2) %>% dplyr::rename(cluster_1 = V1,
                                               cluster_2 = V2)

sites_trout <- sites %>% filter(site %in% unique(k2_trout[[1]]))

k2_trout_structure <- k2_trout %>% 
  dplyr::mutate(id = as.integer(row_number())) %>% 
  gather('pop', 'prob', cluster_1:cluster_2) %>% 
  group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k2_trout_structure$site <- ordered(k2_trout_structure$site,
                                   levels = rev(c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","BRO")))


#-> K=10

k10_trout <- read.table("data/ngs_admix/trout/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_10_2.qopt")
k10_trout$ind <- ind_trout$V1

k10_trout %<>% left_join(trout_apclus,by=c("ind"="id"))

k10_trout %<>% mutate(site = case_when(str_detect(ind,"ROG") ~ "ROG",
                                      str_detect(ind,"CHIm") ~ "KAA",
                                      str_detect(ind,"BIA") ~ "BIA",
                                      str_detect(ind,"LAG") ~ "LAG",
                                      str_detect(ind,"AQU") ~ "AQU",
                                      str_detect(ind,"BEA") ~ "BEA",
                                      str_detect(ind,"SAB") ~ "SAB",
                                      str_detect(ind,"MOA") ~ "SAB",
                                      str_detect(ind,"REN") ~ "REN",
                                      str_detect(ind,"CON") ~ "CON",
                                      str_detect(ind,"EAS") ~ "EAS",
                                      str_detect(ind,"FIS") ~ "EAS",
                                      str_detect(ind,"WASm") & cluster == "BRO" ~ "BRO",
                                      str_detect(ind,"WASm") & cluster == "EAS" ~ "EAS",
                                      str_detect(ind,"WASm") & cluster == "CON" ~ "CON",
                                      str_detect(ind,"WASm") & cluster == "REN" ~ "REN",
                                      str_detect(ind,"RUP") ~ "RUP",
                                      str_detect(ind,"BRO") ~ "BRO",
                                      str_detect(ind,"NOT") ~ "NOT"))


k10_trout %<>% dplyr::select(site,ind,V1:V10) %>% dplyr::rename(cluster_1 = V1,
                                                  cluster_2 = V2,
                                                  cluster_3 = V3,
                                                  cluster_4 = V4,
                                                  cluster_5 = V5,
                                                  cluster_6 = V6,
                                                  cluster_7 = V7,
                                                  cluster_8 = V8,
                                                  cluster_9 = V9,
                                                  cluster_10 = V10)

sites_trout <- sites %>% filter(site %in% unique(k10_trout[[1]]))


k10_trout_structure <- k10_trout %>% 
  dplyr::mutate(id = as.integer(row_number())) %>% 
  gather('pop', 'prob', cluster_1:cluster_10) %>% 
  group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k10_trout_structure$site <- ordered(k10_trout_structure$site,
                                   levels = rev(c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","FIS","EAS","BRO")))

#------------------------------#
#---> Make Structure plots <---#
#------------------------------#


pal_custom1 <- c("#313695","#a50026")


#---> Structure plot K=2
struct_plot_k2_trout <- ggplot(k2_trout_structure, aes(id, prob, fill = pop)) +
  geom_col(width = 2) +
  scale_fill_manual(values = pal_custom1)+
  facet_grid(~site,scales = "free_x")+
  ylab("Ancestry")+
  xlab("")+
  ggtitle("K=2")+
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
        strip.text = element_text(size=18),
        legend.position = "none")


#---> Structure plot K=10

pal_custom2 <- c("#abdda4", #Beaver
                 "#F67E4BFF", #Eastmain
                 "#4575b4", #LaGrande
                 "#d53e4f", #Broadback
                 "#74add1", #Bia
                 "#FEDA8BFF", #Renoyer
                 "#66c2a5", #Aquatuc
                 "#FDB366FF", #Conn
                 "#313695", #Roggan/Kaapsaoui
                 "#EAECCCFF") #Sabacunica

struct_plot_k10_trout <- ggplot(k10_trout_structure, aes(id, prob, fill = pop)) +
  geom_col(width = 2) +
  scale_fill_manual(values = pal_custom2)+
  facet_grid(~site,scales = "free_x")+
  ylab("Ancestry")+
  xlab("")+
  ggtitle("K=10")+
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
        strip.text = element_text(size=18),
        legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_blank())



#------------------------------#
#---> Make mapmixture plot <---#
#------------------------------#

# Parameters
crs <- 4326
boundary <- c(xmin=-81.75, xmax=-78, ymin=51.1, ymax=54.75) |> transform_bbox(bbox = _, crs)

# Read in world coastlines and transform to CRS
file <- system.file("extdata", "world.gpkg", package = "mapmixture")
world <- st_read(file, quiet = TRUE) |> st_transform(x = _, crs = crs)

# Run mapmixture helper functions to prepare admixture and coordinates data
admixture_df_trout <- standardise_data(k10_trout, type = "admixture") |> transform_admix_data(data = _)
coords_df_trout <- standardise_data(sites_trout, type = "coordinates")
admix_coords_trout <- merge_coords_data(coords_df_trout, admixture_df_trout ) |> transform_df_coords(df = _, crs = crs)


trout_admixmap_k10 <-ggplot() +
  geom_sf(data = world, fill = "grey87",color = "black")+
  geom_sf(data =rog,color = "white",linewidth = 0.3)+
  geom_sf(data =kaa,color="white",linewidth = 0.3)+
  geom_sf(data =bia,color="white",linewidth = 0.3)+
  geom_sf(data =lagr,color="white",linewidth = 0.3)+
  geom_sf(data =aqu,color="white",linewidth = 0.3)+
  geom_sf(data =bea,color="white",linewidth = 0.3)+
  geom_sf(data =sab,color="white",linewidth = 0.3)+
  geom_sf(data =sabl,color="white",linewidth = 0.3)+
  geom_sf(data =ren,color="white",linewidth = 0.3)+
  geom_sf(data =conn,color="white",linewidth = 0.3)+
  geom_sf(data =eas,color="white",linewidth = 0.3)+
  geom_sf(data =jac,color="white",linewidth = 0.3)+
  geom_sf(data =rup,color="white",linewidth = 0.3)+
  geom_sf(data =bro,color="white",linewidth = 0.3)+
  geom_sf(data =not,color="white",linewidth = 0.3)+
  coord_sf(
    xlim = c(boundary[["xmin"]], boundary[["xmax"]]),
    ylim = c(boundary[["ymin"]], boundary[["ymax"]]))+
  add_pie_charts(admix_coords_trout,
                 admix_columns = 4:ncol(admix_coords_trout),
                 lat_column = "lat",
                 lon_column = "lon",
                 pie_colours = pal_custom2,
                 border = 0.4, pie_size = 0.3, opacity = 1)+
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white"),text_cex = 1) +
  north(location = "topleft",
        x.min = boundary[["xmin"]], 
        x.max =boundary[["xmax"]], 
        y.min = boundary[["ymin"]], 
        y.max =boundary[["ymax"]])+
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_text(size = 14), #,angle = 45,vjust = .9,hjust=1),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16))

ggplot() +
  geom_sf(data = world, fill = "grey87",color = "black")+
  geom_sf(data =rog,color = "white",linewidth = 0.3)+
  geom_sf(data =kaa,color="white",linewidth = 0.3)+
  geom_sf(data =bia,color="white",linewidth = 0.3)+
  geom_sf(data =lagr,color="white",linewidth = 0.3)+
  geom_sf(data =aqu,color="white",linewidth = 0.3)+
  geom_sf(data =bea,color="white",linewidth = 0.3)+
  geom_sf(data =sab,color="white",linewidth = 0.3)+
  geom_sf(data =sabl,color="white",linewidth = 0.3)+
  geom_sf(data =ren,color="white",linewidth = 0.3)+
  geom_sf(data =conn,color="white",linewidth = 0.3)+
  geom_sf(data =eas,color="white",linewidth = 0.3)+
  geom_sf(data =jac,color="white",linewidth = 0.3)+
  geom_sf(data =rup,color="white",linewidth = 0.3)+
  geom_sf(data =bro,color="white",linewidth = 0.3)+
  geom_sf(data =not,color="white",linewidth = 0.3)+
  coord_sf(
    xlim = c(boundary[["xmin"]], boundary[["xmax"]]),
    ylim = c(boundary[["ymin"]], boundary[["ymax"]]))+
  add_pie_charts(admix_coords_trout,
                 admix_columns = 4:ncol(admix_coords_trout),
                 lat_column = "lat",
                 lon_column = "lon",
                 pie_colours = pal_custom2,
                 border = 0.6, pie_size = 0.3, opacity = 1)+
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white"),text_cex = 1) +
  north(location = "topleft",
        x.min = boundary[["xmin"]], 
        x.max =boundary[["xmax"]], 
        y.min = boundary[["ymin"]], 
        y.max =boundary[["ymax"]])+
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_text(size = 18,angle = 45,vjust = .9,hjust=1),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 16))

ggsave("figures/mapmixture_trout.pdf",height = 8, width = 6)

#-------------------#
#---> Whitefish <---#
#-------------------#

whitefish_apclus <- read_rds("data/ap_cluster/ap_cluster_whitefish.Rds")

whitefish_apclus %<>% mutate(cluster = case_when(L1 == 1 ~ "BEA",
                                                 L1 == 2 ~ "LAG",
                                                 L1 == 3 ~ "KAA",
                                                 L1 == 4 ~ "EAS-C",
                                                 L1 == 5 ~ "RUP",
                                                 L1 == 6 ~ "NOT",
                                                 L1 == 7 ~ "SAB",
                                                 L1 == 8 ~ "ROG"))

whitefish_apclus %<>% dplyr::select(id,cluster)

#---> Start by making dataframes for the relevant structure plot 

#-> K=2

k2_whitefish <- read.table("data/ngs_admix/whitefish/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_2_1.qopt")

ind_whitefish2 <- read.table("data/ngs_admix/whitefish/ind2.txt")

k2_whitefish$ind <- ind_whitefish2$V1

k2_whitefish %<>% left_join(whitefish_apclus,by=c("ind"="id"))

k2_whitefish %<>% mutate(site = case_when(str_detect(ind,"ROG") ~ "ROG",
                                      str_detect(ind,"CHIm") ~ "KAA",
                                      str_detect(ind,"BIA") ~ "BIA",
                                      str_detect(ind,"LAG") & cluster == "LAG" ~ "LAG",
                                      str_detect(ind,"LAG") & cluster == "BEA" ~ "LAG",
                                      str_detect(ind,"LAG") & cluster == "ROG" ~ "ROG",
                                      str_detect(ind,"AQU") ~ "AQU",
                                      str_detect(ind,"BEA") ~ "BEA",
                                      str_detect(ind,"SAB") ~ "SAB",
                                      str_detect(ind,"MOA") ~ "SAB",
                                      str_detect(ind,"REN") ~ "REN",
                                      str_detect(ind,"TIL") ~ "REN",
                                      str_detect(ind,"CON") ~ "CON",
                                      str_detect(ind,"EAS") ~ "EAS",
                                      str_detect(ind,"FIS") ~ "EAS",
                                      str_detect(ind,"JAC") ~ "JAC",
                                      str_detect(ind,"RUP") ~ "RUP",
                                      str_detect(ind,"BRO") ~ "BRO",
                                      str_detect(ind,"NOT") ~ "NOT"))


k2_whitefish %<>% dplyr::select(site,ind,V1,V2) %>% dplyr::rename(cluster_1 = V1,
                                                cluster_2 = V2)

sites_whitefish <- sites %>% filter(site %in% unique(k2_whitefish[[1]]))

k2_whitefish_structure <- k2_whitefish %>% 
  dplyr::mutate(id = as.integer(row_number())) %>% 
  gather('pop', 'prob', cluster_1:cluster_2) %>% 
  group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k2_whitefish_structure$site <- ordered(k2_whitefish_structure$site,
                                       levels = rev(c("ROG","KAA","BIA","LAG","BEA","SAB","REN","CON","EAS","JAC","RUP","NOT")))


#-> K=8

k8_whitefish <- read.table("data/ngs_admix/whitefish/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_8_18.qopt")

ind_whitefish_2 <- read.table("data/ngs_admix/whitefish/ind2.txt")

k8_whitefish$ind <- ind_whitefish_2$V1

k8_whitefish %<>% left_join(whitefish_apclus,by=c("ind"="id"))

k8_whitefish %<>% mutate(site = case_when(str_detect(ind,"ROG") ~ "ROG",
                                          str_detect(ind,"CHIm") ~ "KAA",
                                          str_detect(ind,"BIA") ~ "BIA",
                                          str_detect(ind,"LAG") & cluster == "LAG" ~ "LAG",
                                          str_detect(ind,"LAG") & cluster == "BEA" ~ "LAG",
                                          str_detect(ind,"LAG") & cluster == "ROG" ~ "ROG",
                                          str_detect(ind,"BEA") ~ "BEA",
                                          str_detect(ind,"SAB") ~ "SAB",
                                          str_detect(ind,"MOA") ~ "SAB",
                                          str_detect(ind,"REN") ~ "REN",
                                          str_detect(ind,"TIL") ~ "REN",
                                          str_detect(ind,"CON") ~ "CON",
                                          str_detect(ind,"EAS") ~ "EAS",
                                          str_detect(ind,"FIS") ~ "EAS",
                                          str_detect(ind,"JAC") ~ "JAC",
                                          str_detect(ind,"RUP") ~ "RUP",
                                          str_detect(ind,"NOT") ~ "NOT"))


k8_whitefish %<>% dplyr::select(site,ind,V1:V8) %>% dplyr::rename(cluster_1 = V1,
                                                      cluster_2 = V2,
                                                      cluster_3 = V3,
                                                      cluster_4 = V4,
                                                      cluster_5 = V5,
                                                      cluster_6 = V6,
                                                      cluster_7 = V7,
                                                      cluster_8 = V8)

sites_whitefish <- sites %>% filter(site %in% unique(k8_whitefish[[1]]))


k8_whitefish_structure <- k8_whitefish %>% 
  dplyr::mutate(id = as.integer(row_number())) %>% 
  gather('pop', 'prob', cluster_1:cluster_8) %>% 
  group_by(id) %>% 
  dplyr::mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k8_whitefish_structure$site <- ordered(k8_whitefish_structure$site,
                                    levels = rev(c("ROG","KAA","BIA","LAG","BEA","SAB","REN","CON","EAS","JAC","RUP","NOT")))


#------------------------------#
#---> Make Structure plots <---#
#------------------------------#

#Structure plot K=2
struct_plot_k2_whitefish <- ggplot(k2_whitefish_structure, aes(id, prob, fill = pop)) +
  geom_col(width = 2) +
  scale_fill_manual(values = pal_custom1)+
  facet_grid(~site,scales = "free_x")+
  ylab("Ancestry")+
  xlab("")+
  ggtitle("K=2")+
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
        strip.text = element_text(size=18),
        legend.position = "none")

#Structure plot K=8
struct_plot_k8_whitefish <- ggplot(k8_whitefish_structure, aes(id, prob, fill = pop)) +
  geom_col(width = 2) +
  scale_fill_manual(values = c("#abdda4", #BEAVER
                               "#f46d43", #RUPERT
                               "#364B9AFF", #ROGGAN
                               "#a50026", #NOTTAWAY
                               "#4575b4", #KAAPSAOUI
                               "#EAECCCFF", #SABACUNICA
                               "#74add1", #LAGRANDE
                               "#fdae61"))+
  facet_grid(~site,scales = "free_x")+
  ylab("Ancestry")+
  xlab("")+
  ggtitle("K=8")+
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
        strip.text = element_text(size=18),
        legend.position = "none",strip.background = element_blank(),
        strip.text.x = element_blank())


#------------------------------#
#---> Make mapmixture plot <---#
#------------------------------#

# Parameters
crs <- 4326
boundary <- c(xmin=-81.75, xmax=-78, ymin=51.1, ymax=54.75) |> transform_bbox(bbox = _, crs)

# Read in world coastlines and transform to CRS
file <- system.file("extdata", "world.gpkg", package = "mapmixture")
world <- st_read(file, quiet = TRUE) |> st_transform(x = _, crs = crs)

# Run mapmixture helper functions to prepare admixture and coordinates data
admixture_df_wf <- standardise_data(k8_whitefish, type = "admixture") |> transform_admix_data(data = _)
coords_df_wf <- standardise_data(sites_whitefish, type = "coordinates")
admix_coords_wf <- merge_coords_data(coords_df_wf, admixture_df_wf ) |> transform_df_coords(df = _, crs = crs)


ggplot() +
  geom_sf(data = world, fill = "grey87",color = "black")+
  geom_sf(data =rog,color = "white",linewidth = 0.3)+
  geom_sf(data =kaa,color="white",linewidth = 0.3)+
  geom_sf(data =bia,color="white",linewidth = 0.3)+
  geom_sf(data =lagr,color="white",linewidth = 0.3)+
  geom_sf(data =bea,color="white",linewidth = 0.3)+
  geom_sf(data =sab,color="white",linewidth = 0.3)+
  geom_sf(data =sabl,color="white",linewidth = 0.3)+
  geom_sf(data =ren,color="white",linewidth = 0.3)+
  geom_sf(data =conn,color="white",linewidth = 0.3)+
  geom_sf(data =eas,color="white",linewidth = 0.3)+
  geom_sf(data =jac,color="white",linewidth = 0.3)+
  geom_sf(data =rup,color="white",linewidth = 0.3)+
  geom_sf(data =bro,color="white",linewidth = 0.3)+
  geom_sf(data =not,color="white",linewidth = 0.3)+
  coord_sf(
    xlim = c(boundary[["xmin"]], boundary[["xmax"]]),
    ylim = c(boundary[["ymin"]], boundary[["ymax"]]))+
  add_pie_charts(admix_coords_wf,
                 admix_columns = 4:ncol(admix_coords_wf),
                 lat_column = "lat",
                 lon_column = "lon",
                 pie_colours = c("#abdda4", #BEAVER
                                 "#f46d43", #RUPERT
                                 "#364B9AFF", #ROGGAN
                                 "#a50026", #NOTTAWAY
                                 "#4575b4", #KAAPSAOUI
                                 "#EAECCCFF", #SABACUNICA
                                 "#74add1", #LAGRANDE
                                 "#fdae61"),
                 border = 0.6, pie_size = 0.3, opacity = 1)+
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white"),text_cex = 1) +
  north(location = "topleft",
        x.min = boundary[["xmin"]], 
        x.max =boundary[["xmax"]], 
        y.min = boundary[["ymin"]], 
        y.max =boundary[["ymax"]])+
  theme(
    legend.title = element_blank(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(colour = "black", fill = NA),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 16))

ggsave("figures/mapmixture_whitefish.pdf",height = 8, width = 6)


#--------------------#
#---> FINAL PLOT <---#
#--------------------#

structure_final <- (struct_plot_k2_whitefish/struct_plot_k8_whitefish/struct_plot_k2_trout/struct_plot_k10_trout)

ggsave(plot =structure_final,"figures/structure_whitefish_trout.png",width = 12)


#---> Trash

#Trout

tmp <- read.table("data/ngs_admix/trout/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_10_22.qopt")
tmp$ind <- ind_trout$V1

tmp %<>% mutate(site = case_when(str_detect(ind,"ROG") ~ "ROG",
                                       str_detect(ind,"CHIm") ~ "KAA",
                                       str_detect(ind,"BIA") ~ "BIA",
                                       str_detect(ind,"LAG") ~ "LAG",
                                       str_detect(ind,"AQU") ~ "AQU",
                                       str_detect(ind,"BEA") ~ "BEA",
                                       str_detect(ind,"SAB") ~ "SAB",
                                       str_detect(ind,"MOA") ~ "SAB",
                                       str_detect(ind,"REN") ~ "REN",
                                       str_detect(ind,"CON") ~ "CON",
                                       str_detect(ind,"EAS") ~ "EAS",
                                       str_detect(ind,"FIS") ~ "EAS",
                                       str_detect(ind,"WASm") ~ "BRO",
                                       str_detect(ind,"RUP") ~ "RUP",
                                       str_detect(ind,"BRO") ~ "BRO",
                                       str_detect(ind,"NOT") ~ "NOT")) %>% 
  dplyr::select(ind,site)

write.table(tmp,"data/ngs_admix/trout/wgs_assign_pop_id.txt",quote=F,row.names=F,col.names = F,sep="\t")

#Whitefish

tmp2 <- read.table("data/ngs_admix/whitefish/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_8_2.qopt")
tmp2$ind <- ind_whitefish$V1

tmp2 %<>% mutate(site = case_when(str_detect(ind,"ROG") ~ "ROG",
                                 str_detect(ind,"CHIm") ~ "KAA",
                                 str_detect(ind,"BIA") ~ "BIA",
                                 str_detect(ind,"LAG") ~ "LAG",
                                 str_detect(ind,"BEA") ~ "BEA",
                                 str_detect(ind,"SAB") ~ "SAB",
                                 str_detect(ind,"MOA") ~ "SAB",
                                 str_detect(ind,"REN") ~ "REN",
                                 str_detect(ind,"TIL") ~ "REN",
                                 str_detect(ind,"CON") ~ "CON",
                                 str_detect(ind,"EAS") ~ "EAS",
                                 str_detect(ind,"FIS") ~ "EAS",
                                 str_detect(ind,"JAC") ~ "JAC",
                                 str_detect(ind,"RUP") ~ "RUP",
                                 str_detect(ind,"NOT") ~ "NOT")) %>% 
  dplyr::select(ind,site)

write.table(tmp2,"data/ngs_admix/whitefish/wgs_assign_pop_id.txt",quote=F,row.names=F,col.names = F,sep="\t")





#TRASH


lag_tmp <- k8_whitefish %>% filter(site== "LAG") %>% select(ind)

rog_tmp <- k8_whitefish %>% filter(site == "ROG") %>% select(ind)

write.table(lag_tmp,"~/Desktop/lag_tmp.tsv",quote = F,row.names = F,col.names = F)

write.table(rog_tmp,"~/Desktop/rog_tmp.tsv",quote = F,row.names = F,col.names = F)


