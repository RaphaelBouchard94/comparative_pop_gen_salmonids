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

#---> Whitefish

#Select individual to equalize effective sample size per cluster

ne_whitefish <- read.table("data/wgs_assign/whitefish/cocl_456_ind_230697_pruned_singleton_snp.ne_ind.txt",header=F)

ind_whitefish <- read.table("data/wgs_assign/whitefish/whitefish_wgs_assign_pop_id.txt")

colnames(ind_whitefish) <- c("id","cluster")

ind_whitefish$ne_whitefish <- ne_whitefish$V1

ind_whitefish %>% group_by(cluster) %>% dplyr::summarize(sum(ne_whitefish))

#Format admixture file to get least admixed individuals for source pop

k8_whitefish <- read.table("data/ngs_admix/whitefish/all_maf0.05_pctind0.8_maxdepth10_pruned_singletons_8_2.qopt")

k8_whitefish$id <- ind_whitefish$id

k8_whitefish %<>% left_join(ind_whitefish %>% dplyr::select(id,cluster), by = "id")

k8_whitefish %<>% dplyr::select(cluster,id,V1:V8) 

k8_whitefish$best_k <- colnames(k8_whitefish)[apply(k8_whitefish,1,which.max)]

k8_whitefish %<>% 
  rowwise() %>% dplyr::mutate(best_k_admix = max(c_across(starts_with("V")), na.rm = TRUE))

k8_whitefish %<>% dplyr::select(id,best_k,best_k_admix)

#Merge admixture and effective sample size 

ind_whitefish %<>% left_join(k8_whitefish,by="id")

ind_whitefish_select <- ind_whitefish %>% group_by(cluster) %>% arrange(desc(best_k_admix))


# Function to find the desired number of observations
select_observations <- function(df_group, constant) {
  selected_obs <- df_group %>%
    mutate(cumsum_value = cumsum(ne_whitefish)) %>%
    filter(cumsum_value <= constant) %>%
    slice(1:n())  # Select all rows until the condition is met
  return(selected_obs)
}

essbp_samples <- ind_whitefish_select %>% 
  group_split(cluster) %>% 
  map_df(~ select_observations(.x, 14))

essbp_samples %>% group_by(cluster) %>% dplyr::count()

essbp_samples %>% group_by(cluster) %>% dplyr::summarise(sum_ne = sum(ne_whitefish))

id_essbp_orderedbam <- ind_whitefish %>% filter(id %in% essbp_samples$id)

#write.table(id_essbp_orderedbam %>% dplyr::select(id),"data/wgs_assign/whitefish/eesbp_samples.txt",quote=F,col.names = F,row.names = F)

#write.table(id_essbp_orderedbam %>% dplyr::select(id,cluster),"data/wgs_assign/whitefish/eesbp_samples_cluster.txt",quote=F,col.names = F,row.names = F,sep="\t")

essbp <- id_essbp_orderedbam %>% dplyr::select(id,true_origin)

#--------------#

loo_whitefish <- read.table("data/wgs_assign/whitefish/cocl_98_ind_eesbp_230697_pruned_singleton_snp.pop_like_LOO.txt")
pop_order_whitefish <- read.table("data/wgs_assign/whitefish/cocl_98_ind_eesbp_230697_pruned_singleton_snp.pop_names.txt")

colnames(loo_whitefish) <- pop_order_whitefish$V1

#Estimate origin of essbp dataset individuals using LOD

loo_whitefish %<>%
  rowwise() %>%
  dplyr::mutate(origin = paste0(names(.)[c_across() == max(c_across())], collapse = '_'))

#Add assignment info to sampling site

id_essbp_orderedbam$true_origin <- loo_whitefish$origin

id_essbp_orderedbam %<>% dplyr::mutate(migrant = if_else(cluster == true_origin, "no","yes"))

id_essbp_orderedbam %>% group_by(migrant) %>% dplyr::count()

#4 inconsistencies (between BEA-SAB and NOT-RUP)

#Assignment accuracy with the ESSBP dataset 96%

#---> Estimate assignment of remaining of the dataset

#write.table(ind_whitefish %>% filter(!id %in% essbp_samples$id) %>% pull(id), "data/wgs_assign/whitefish/non_essbp_samples.txt",quote=F,col.names = F,row.names = F)

likehood_trout_non_essbp <- read.table("data/wgs_assign/whitefish/cocl_358_remaining_samples_assignment.pop_like.txt")

colnames(likehood_trout_non_essbp) <- pop_order_whitefish$V1

#Estimate origin of non essbp dataset individuals using LOD

likehood_trout_non_essbp %<>%
  rowwise() %>%
  dplyr::mutate(true_origin = paste0(names(.)[c_across() == max(c_across())], collapse = '_'))

#add id and name

likehood_trout_non_essbp$id <- ind_whitefish %>% filter(!id %in% essbp_samples$id) %>% pull(id)
likehood_trout_non_essbp$origin <- ind_whitefish %>% filter(!id %in% essbp_samples$id) %>% pull(cluster)

#Test if assignment matches

likehood_trout_non_essbp %<>% dplyr::mutate(migrant = if_else(true_origin == origin, "no","yes"))

likehood_trout_non_essbp %>% group_by(migrant) %>% dplyr::count()

#89 inconsistencies

non_essbp <- likehood_trout_non_essbp %>% dplyr::select(id,true_origin)

#-------------------------------------#
#---> Estimate number of migrants <---#
#-------------------------------------#

all_samples <- rbind(essbp,non_essbp)

all_samples %<>% dplyr::mutate(site = case_when(str_detect(id,"ROG") ~ "Roggan",
                                                str_detect(id,"CHIm") ~ "Kaapsaoui",
                                                str_detect(id,"BIA") ~ "La_Grande",
                                                str_detect(id,"LAG") ~ "La_Grande",
                                                str_detect(id,"BEA") ~ "Beaver",
                                                str_detect(id,"SAB") ~ "Sabacunica",
                                                str_detect(id,"MOA") ~ "Sabacunica",
                                                str_detect(id,"REN") ~ "Eastmain",
                                                str_detect(id,"TIL") ~ "Eastmain",
                                                str_detect(id,"CON") ~ "Eastmain",
                                                str_detect(id,"EAS") ~ "Eastmain",
                                                str_detect(id,"FIS") ~ "Eastmain",
                                                str_detect(id,"JAC") ~ "Eastmain",
                                                str_detect(id,"RUP") ~ "Rupert",
                                                str_detect(id,"NOT") ~ "Nottaway"))

all_samples %<>% dplyr::mutate(migrant = if_else(site == true_origin, "no","yes"))

all_samples %>% dplyr::group_by(migrant) %>% dplyr::count()

#172 dispersal event

all_samples %>% dplyr::group_by(true_origin) %>% dplyr::count()










