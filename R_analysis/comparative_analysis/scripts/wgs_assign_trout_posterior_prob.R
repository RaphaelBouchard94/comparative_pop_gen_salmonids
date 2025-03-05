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

pop_order_trout <- read.table("data/wgs_assign/trout/safo_254_essbp_ind_118861_pruned_singleton_snp.pop_names.txt")

#---> LOO

ind_id <- read.table("data/wgs_assign/trout/eesbp_samples_cluster.txt")

#List les fichiers LOO

files <- list.files(path = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/wgs_assign/trout/posterior_prob/", 
                    pattern = "LOO",
                    full.names = T)

#Importe les fichiers dans une liste puis fusionne la liste

file_list = list()
for (i in 1:length(files)){
  a <-  read.table(files[i])
  file_list[[i]] = a
}

get_best_lik <- function(setx){
  #Add pop names
  colnames(setx) <- pop_order_trout$V1
  
  #Estimate origin of essbp dataset individuals using LOD
  setx %<>%
    rowwise() %>%
    dplyr::mutate(origin = paste0(names(.)[c_across() == max(c_across())], collapse = '_'))
  
}

#Get best lik for each element of the list
file_list_best_lik <- lapply(file_list,get_best_lik)

#Merge the column origin of all element of the list into a single data.frame
post_prob_essbp <- as.data.frame(lapply(colnames(file_list_best_lik[[1]])[11], function(x) sapply(file_list_best_lik, `[[`, x))[[1]])

#add id
post_prob_essbp$id <- ind_id$V1

#Assign individual based on posterior probability
post_prob_essbp <- pivot_longer(post_prob_essbp,cols = V1:V10,names_to = "rep",values_to = "best_lik")

post_prob_essbp %<>% group_by(id) %>% dplyr::count(best_lik) %>% transmute(best_lik,post_prob = n/10) 

post_prob_essbp %<>% mutate(assigned = if_else(post_prob >= 0.8,"yes","no"))

post_prob_essbp %<>% group_by(id) %>% slice(which.max(post_prob))

post_prob_essbp %>% group_by(assigned) %>% dplyr::count()

post_prob_essbp %>% group_by(assigned,best_lik) %>% dplyr::count()


#---> Basic assignment

#get id
non_essbp <- read.table("data/wgs_assign/trout/non_essbp_samples.txt")

#get files
files_assign <- list.files(path = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/wgs_assign/trout/posterior_prob/", 
                           pattern = "pop_like.txt",
                           full.names = T)

#Importe les fichiers dans une liste puis fusionne la liste

file_list2 = list()
for (i in 1:length(files_assign)){
  a <-  read.table(files_assign[i])
  file_list2[[i]] = a
}

#Get best lik for each element of the list
file_list_best_lik_assign <- lapply(file_list2,get_best_lik)

#Merge the column origin of all element of the list into a single data.frame
post_prob_non_essbp <- as.data.frame(lapply(colnames(file_list_best_lik_assign[[1]])[11], function(x) sapply(file_list_best_lik_assign, `[[`, x))[[1]])

#add id
post_prob_non_essbp$id <- non_essbp$V1

#Assign individual based on posterior probability
post_prob_non_essbp <- pivot_longer(post_prob_non_essbp,cols = V1:V10,names_to = "rep",values_to = "best_lik")

post_prob_non_essbp %<>% group_by(id) %>% dplyr::count(best_lik) %>% transmute(best_lik,post_prob = n/10) 

post_prob_non_essbp %<>% mutate(assigned = if_else(post_prob >= 0.8,"yes","no"))

post_prob_non_essbp %<>% group_by(id) %>% slice(which.max(post_prob))

post_prob_non_essbp %>% group_by(assigned) %>% dplyr::count()

#--> Merge essbp and non essbp in a single file

trout_assignment <- rbind(post_prob_essbp,post_prob_non_essbp)

write.table(trout_assignment,"data/dispersal_model/trout_assignment.txt",quote = F,row.names = F)








