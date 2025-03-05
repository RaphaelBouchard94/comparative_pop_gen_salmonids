rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(patchwork)
library(plyr)

#---------------#
#----- Data ----#
#---------------#

#---> TROUT

#Select individual to equalize effective sample size per cluster

ne_trout <- read.table("data/wgs_assign/trout/safo_362_ind_118861_pruned_singleton_snp.ne_ind.txt",header=F)

ind_trout <- read.table("data/wgs_assign/trout/wgs_assign_pop_id.txt")

colnames(ind_trout) <- c("id","cluster")

ind_trout$ne_trout <- ne_trout$V1

#Format admixture file to get least admixed individuals for source pop

k10_trout <- read.table("data/ngs_admix/trout/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_10_1.qopt")

k10_trout$id <- ind_trout$id

k10_trout %<>% left_join(ind_trout %>% dplyr::select(id,cluster), by = "id")

k10_trout %<>% dplyr::select(cluster,id,V1:V10) 

k10_trout$best_k <- colnames(k10_trout)[apply(k10_trout,1,which.max)]

k10_trout %<>% 
  rowwise() %>% dplyr::mutate(best_k_admix = max(c_across(starts_with("V")), na.rm = TRUE))

k10_trout %<>% dplyr::select(id,best_k,best_k_admix)

#Merge admixture and effective sample size 

ind_trout %<>% left_join(k10_trout,by="id")

ind_trout_select <- ind_trout %>% group_by(cluster) %>% arrange(desc(best_k_admix))


# Function to find the desired number of observations
select_observations <- function(df_group, constant) {
  selected_obs <- df_group %>%
    mutate(cumsum_value = cumsum(ne_trout)) %>%
    filter(cumsum_value <= constant) %>%
    slice(1:n())  # Select all rows until the condition is met
  return(selected_obs)
}

essbp_samples <- ind_trout_select %>% 
  group_split(cluster) %>% 
  map_df(~ select_observations(.x, 15))

essbp_samples %>% group_by(cluster) %>% dplyr::count()

 essbp_samples %>% group_by(cluster) %>% dplyr::summarise(sum_ne = sum(ne_trout))

id_essbp_orderedbam <- ind_trout %>% filter(id %in% essbp_samples$id)

#write.table(id_essbp_orderedbam %>% dplyr::select(id),"data/wgs_assign/trout/eesbp_samples.txt",quote=F,col.names = F,row.names = F)

#write.table(id_essbp_orderedbam %>% dplyr::select(id,cluster),"data/wgs_assign/trout/eesbp_samples_cluster.txt",quote=F,col.names = F,row.names = F,sep="\t")


#--------------#

loo_trout <- read.table("data/wgs_assign/trout/safo_254_essbp_ind_118861_pruned_singleton_snp.pop_like_LOO.txt")
pop_order_trout <- read.table("data/wgs_assign/trout/safo_254_essbp_ind_118861_pruned_singleton_snp.pop_names.txt")

colnames(loo_trout) <- pop_order_trout$V1

#Estimate origin of essbp dataset individuals using LOD

loo_trout %<>%
  rowwise() %>%
  dplyr::mutate(origin = paste0(names(.)[c_across() == max(c_across())], collapse = '_'))

#Add assignment info to sampling site

id_essbp_orderedbam$true_origin <- loo_trout$origin

id_essbp_orderedbam %<>% dplyr::mutate(migrant = if_else(cluster == true_origin, "no","yes"))

id_essbp_orderedbam %>% group_by(migrant) %>% dplyr::count()

essbp <- id_essbp_orderedbam %>% dplyr::select(id,true_origin)

#6 errors


#---> Estimate assignment of remaining of the dataset

#write.table(ind_trout %>% filter(!id %in% essbp_samples$id) %>% pull(id), "data/wgs_assign/trout/non_essbp_samples.txt",quote=F,col.names = F,row.names = F)

likehood_trout_non_essbp <- read.table("data/wgs_assign/trout/safo_120_remaining_samples_assignment.pop_like.txt")

colnames(likehood_trout_non_essbp) <- pop_order_trout$V1

#Estimate origin of non essbp dataset individuals using LOD

likehood_trout_non_essbp %<>%
  rowwise() %>%
  dplyr::mutate(true_origin = paste0(names(.)[c_across() == max(c_across())], collapse = '_'))

#add id and name

likehood_trout_non_essbp$id <- ind_trout %>% filter(!id %in% essbp_samples$id) %>% pull(id)
likehood_trout_non_essbp$origin <- ind_trout %>% filter(!id %in% essbp_samples$id) %>% pull(cluster)

#Test if assignment matches

likehood_trout_non_essbp %<>% dplyr::mutate(migrant = if_else(true_origin == origin, "no","yes"))

likehood_trout_non_essbp %>% group_by(migrant) %>% dplyr::count()

#21 inconsistencies with APcluster + NGSadmix assignment

non_essbp <- likehood_trout_non_essbp %>% dplyr::select(id,true_origin)

#Estimate number of migrants

all_samples <- rbind(essbp,non_essbp)


all_samples %<>% dplyr::mutate(site = case_when(str_detect(id,"ROG") ~ "Kaapsaoui",
                                               str_detect(id,"CHIm") ~ "Kaapsaoui",
                                               str_detect(id,"BIA") ~ "Biagadwshi",
                                               str_detect(id,"LAG") ~ "La_Grande",
                                               str_detect(id,"AQU") ~ "Aquatuc",
                                               str_detect(id,"BEA") ~ "Beaver",
                                               str_detect(id,"SAB") ~ "Sabacunica",
                                               str_detect(id,"MOA") ~ "Sabacunica",
                                               str_detect(id,"REN") ~ "Renoyer",
                                               str_detect(id,"CON") ~ "Conn",
                                               str_detect(id,"EAS") ~ "Eastmain",
                                               str_detect(id,"FIS") ~ "Eastmain",
                                               str_detect(id,"WASm") ~ "Broadback",
                                               str_detect(id,"RUP") ~ "Rupert",
                                               str_detect(id,"BRO") ~ "Broadback",
                                               str_detect(id,"NOT") ~ "Nottaway"))




all_samples %<>% dplyr::mutate(migrant = if_else(site == true_origin, "no","yes"))

all_samples %>% dplyr::group_by(migrant) %>% dplyr::count()

#75 dispersal event












