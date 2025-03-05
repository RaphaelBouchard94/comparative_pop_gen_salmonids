rm(list = ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#-----------------#
#---- LIBRARY ----#
#-----------------#

library(tidyverse)
library(magrittr)
library(marmap)
library(DescTools)
library(reshape2)

#-----------------#
#---- ANALYSIS ---#
#-----------------#

#--------------#
#---- TROUT ---#
#--------------#

#---> Get pairwise distance for brook trout (Remove Rupert and Nottaway)

JB <- getNOAA.bathy(lon1 = -82.084, lon2 = -78,
                    lat1 = 51.021, lat2 = 55, resolution = 1)

sites <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/pop_dispersal.txt",header=T)

sites_bt <- plyr::arrange(sites, pop_code, desc(pop_code)) %>% filter(!pop_code %in% c("JAC","RUP","NOT"))

marsites_bt <- sites_bt %>% dplyr::select(long,lat) %>% dplyr::rename(x=long,y=lat)

trans1 <- trans.mat(JB,min.depth = 0)

dist_bt <- lc.dist(trans1,marsites_bt,res="dist")

dist_bt

mat_dist_bt<-as.matrix(dist_bt)

dfdist_bt <- setNames(melt(mat_dist_bt), c('pop1', 'pop2', 'values'))


#---> Import assignment for brook trout

trout_assignment <- read.table("data/dispersal_model/trout_assignment.txt",header=T)

trout_assignment %<>% mutate(site = case_when(str_detect(id,"ROG") ~ "ROG",
                                              str_detect(id,"CHIm") ~ "KAA",
                                              str_detect(id,"BIA") ~ "BIA",
                                              str_detect(id,"LAG") ~ "LAG",
                                              str_detect(id,"AQU") ~ "AQU",
                                              str_detect(id,"BEA") ~ "BEA",
                                              str_detect(id,"SAB") ~ "SAB",
                                              str_detect(id,"MOA") ~ "SAB",
                                              str_detect(id,"REN") ~ "REN",
                                              str_detect(id,"CON") ~ "CON",
                                              str_detect(id,"EAS") ~ "EAS",
                                              str_detect(id,"FIS") ~ "EAS",
                                              str_detect(id,"BRO") ~ "BRO"))


trout_assignment %<>% mutate(best_lik2 = case_when(str_detect(best_lik,"Kaapsaoui") & site == "ROG" ~ "ROG",
                                                   str_detect(best_lik,"Kaapsaoui") & site == "KAA" ~ "KAA",
                                                   str_detect(best_lik,"Kaapsaoui") & site == "BIA" ~ "BIA",
                                                   str_detect(best_lik,"Biagadwshi") ~ "BIA",
                                                   str_detect(best_lik,"La_Grande") ~ "LAG",
                                                   str_detect(best_lik,"Aquatuc") ~ "AQU",
                                                   str_detect(best_lik,"Beaver") ~ "BEA",
                                                   str_detect(best_lik,"Sabacunica") ~ "SAB",
                                                   str_detect(best_lik,"Renoyer") ~ "REN",
                                                   str_detect(best_lik,"Conn") ~ "CON",
                                                   str_detect(best_lik,"Eastmain") ~ "EAS",
                                                   str_detect(best_lik,"Broadback") ~ "BRO"))

trout_assignment %<>% mutate(site = if_else(str_detect(id,"WAS"),best_lik2,site))


#Create a column which identifies individual as migrant

trout_assignment %<>% mutate(migrant = if_else(best_lik2 == site,"no","yes"))

trout_assignment %>% group_by(migrant) %>% dplyr::count(post_prob >= 0.8)

#--->
#Dispersal analysis

Nmax <- trout_assignment %>% filter(migrant=="yes",post_prob >= 0.8) %>% dplyr::count() %>% pull()

#Create vector of probabilities for random draw
probf <- as.vector(table(trout_assignment$site)/sum(table(trout_assignment$site)))

#expected random distribution of geographical distance 
#between source and misassigned rivers was generated 
#by drawing without replacement 10 000 pairs of populations

#rdraw <- lapply(1:10000, function(i) lapply(1:Nmax, function(i) sample(1:11, 2, replace = F, prob = probf)))

#saveRDS(rdraw, file = "data/dispersal_model/simulation_trout.rds")

rdraw <- read_rds("data/dispersal_model/simulation_trout.rds")

# Convert the list to a data frame
rdraw_df <- do.call(rbind, rdraw)

# Convert to a data frame and name the columns
rdraw_df <- as.data.frame(matrix(unlist(rdraw), ncol = 2, byrow = TRUE))
colnames(rdraw_df) <- c("pop1", "pop2")

# Add a column specifying the list number
rdraw_df$list_number <- rep(1:10000, each = Nmax)

# Add distance value

rdraw_dist <- left_join(rdraw_df, dfdist_bt, by = c("pop1","pop2"))

#hist(rdraw_dist$values)

# Add distance class

rdraw_dist %<>% mutate(dist_class = case_when(values < 50 ~ "]0-50[",
                                              values >= 50 & values < 100 ~ "[50-100[",
                                              values >= 100 & values < 150 ~ "[100-150[",
                                              values >= 150 & values < 200 ~ "[150-200[",
                                              values >= 200 & values < 250 ~ "[200-250[",
                                              values >= 250 & values < 300 ~ "[250-300[",
                                              values >= 300 & values <= 355 ~ "[300-355["))

# Count number of individual falling in each class per list

rdraw_dist_count <- rdraw_dist %>% group_by(list_number) %>% dplyr::count(dist_class)

rdraw_dist_count$n <- rdraw_dist_count$n/Nmax

rdraw_dist_count %<>% group_by(dist_class) %>% dplyr::summarise(avg_n = mean(n),
                                                                lower_95ci = quantile(n, probs = 0.025),
                                                                upper_95ci = quantile(n,probs = 0.975))


#Relevel factors
rdraw_dist_count$dist_class <- as.factor(rdraw_dist_count$dist_class)

rdraw_dist_count$dist_class <- factor(rdraw_dist_count$dist_class,
                                      levels = c("]0-50[", "[50-100[", "[100-150[", "[150-200[", "[200-250[","[250-300[","[300-355["))


#Plot null distribution
ggplot(rdraw_dist_count) +
  geom_bar( aes(x=dist_class, y=avg_n), stat="identity", fill="white",color = "black", alpha=0.5) +
  geom_pointrange( aes(x=dist_class, y=avg_n, ymin=lower_95ci, ymax=upper_95ci), colour="darkred", alpha=0.9, size=1.3)+
  theme_bw()

#Average number of dispersers under random dispersal simulation

avg_dist <- rdraw_dist %>% group_by(list_number) %>% dplyr::summarise(mean_dist = mean(values))

hist(avg_dist$mean_dist)

mean(avg_dist$mean_dist)
quantile(avg_dist$mean_dist, probs = 0.025)
quantile(avg_dist$mean_dist, probs = 0.975)


# The average assignment distance of 57 individual brook trout randomly 
# dispersing over the 11 sampled populations by simulations was 
# 127 km [95%CI = (19 – 147)].

#--->
# Add observed distribution

trout_assignment$pop1 <- as.numeric(as.factor(trout_assignment$site))
trout_assignment$pop2 <- as.numeric(as.factor(trout_assignment$best_lik2))

trout_disp <- trout_assignment %>% 
  filter(migrant=="yes",post_prob >= 0.8) %>% 
  dplyr::select(id,pop1,pop2) 


trout_disp %<>% left_join(dfdist_bt)

hist(trout_disp$values)

mean(trout_disp$values)

t.test(log(trout_disp$values),log(avg_dist$mean_dist))
  
#In trout, the 47 misassigned individuals were assigned at 
#an average distance of 37 km from their sampling river

trout_disp  %<>% mutate(dist_class = case_when(values < 50 ~ "]0-50[",
                                           values >= 50 & values < 100 ~ "[50-100[",
                                           values >= 100 & values < 150 ~ "[100-150[",
                                           values >= 150 & values < 200 ~ "[150-200[",
                                           values >= 200 & values < 250 ~ "[200-250[",
                                           values >= 250 & values < 300 ~ "[250-300[",
                                           values >= 300 & values <= 355 ~ "[300-355["))

df_disp_obs <- trout_disp %>% dplyr::count(dist_class)

df_disp_obs$n <- df_disp_obs$n/Nmax

df_disp_obs <- rbind(df_disp_obs,c("[150-200[",0))
df_disp_obs <- rbind(df_disp_obs,c("[200-250[",0))
df_disp_obs <- rbind(df_disp_obs,c("[250-300[",0))
df_disp_obs <- rbind(df_disp_obs,c("[300-355[",0))

df_disp_obs$dist_class <- as.factor(df_disp_obs$dist_class)

df_disp_obs$n <- as.numeric(df_disp_obs$n)

df_disp_obs$dist_class <- factor(df_disp_obs$dist_class,
                                 levels = c("]0-50[", "[50-100[", "[100-150[", "[150-200[", "[200-250[","[250-300[","[300-355["))

ggplot() +
  geom_ribbon(data = rdraw_dist_count,aes(ymin = lower_95ci, ymax = upper_95ci,x= dist_class,group=1), fill = "grey70",alpha = 0.4,colour = "black",linetype=2) +
  geom_line(data = rdraw_dist_count,aes(y = avg_n, x = dist_class,group=1))+
  geom_point(data = df_disp_obs, aes(x = dist_class, y = n),color="darkred")+
  geom_line(data = df_disp_obs, aes(x = dist_class, y = n,group=1),color="darkred")+
  ylab("Proportion of dispersers")+
  xlab("\nDispersal distance (km)")+
  ylim(c(0,1))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle = 45,vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("figures/dispersal_trout.png",width =7, height = 7)


df_disp_obs %<>% left_join(rdraw_dist_count, by = "dist_class")

#GTEST whitefish

gt_disp_obs_bt <- trout_disp %>% dplyr::count(dist_class)

gt_disp_obs_bt <- rbind(gt_disp_obs_bt,c("[150-200[",0))
gt_disp_obs_bt <- rbind(gt_disp_obs_bt,c("[200-250[",0))
gt_disp_obs_bt <- rbind(gt_disp_obs_bt,c("[250-300[",0))
gt_disp_obs_bt <- rbind(gt_disp_obs_bt,c("[300-355[",0))

rdraw_dist_count$n <- rdraw_dist_count$avg_n*47

gt_bt <- gt_disp_obs_bt%>% left_join(rdraw_dist_count, by = "dist_class")

gt_bt_f <- rbind(as.numeric(gt_bt$n.x),gt_bt$n.y)

GTest(gt_bt_f, correct = "none")

#------------------#
#---- WHITEFISH ---#
#------------------#

#---> Get pairwise distance for brook trout (Remove Rupert and Nottaway)

sites_wf <- plyr::arrange(sites, pop_code, desc(pop_code)) %>% filter(!pop_code %in% c("AQU","BRO"))

marsites_wf <- sites_wf %>% dplyr::select(long,lat) %>% dplyr::rename(x=long,y=lat)

dist_wf <- lc.dist(trans1,marsites_wf,res="dist")

dist_wf

mat_dist_wf<-as.matrix(dist_wf)

dfdist_wf <- setNames(melt(mat_dist_wf), c('pop1', 'pop2', 'values'))


#---> Import assignment for whitefish

whitefish_assignment <- read.table("data/dispersal_model/whitefish_assignment.txt",header=T)

whitefish_assignment %<>% mutate(site = case_when(str_detect(id,"ROG") ~ "ROG",
                                                  str_detect(id,"CHIm") ~ "KAA",
                                                  str_detect(id,"BIA") ~ "LAG",
                                                  str_detect(id,"LAG") & best_lik == "La_Grande" ~ "LAG",
                                                  str_detect(id,"LAG") & best_lik == "Roggan" ~ "ROG",
                                                  str_detect(id,"LAG") & best_lik == "Sabacunica" ~ "LAG",
                                                  str_detect(id,"BEA") ~ "BEA",
                                                  str_detect(id,"SAB") ~ "SAB",
                                                  str_detect(id,"MOA") ~ "SAB",
                                                  str_detect(id,"REN") ~ "EAS",
                                                  str_detect(id,"TIL") ~ "EAS",
                                                  str_detect(id,"CON") ~ "EAS",
                                                  str_detect(id,"EAS") ~ "EAS",
                                                  str_detect(id,"FIS") ~ "EAS",
                                                  str_detect(id,"JAC") ~ "EAS",
                                                  str_detect(id,"RUP") ~ "RUP",
                                                  str_detect(id,"NOT") ~ "NOT"))

whitefish_assignment %<>% mutate(best_lik2 = case_when(str_detect(best_lik,"Roggan") ~ "ROG",
                                                       str_detect(best_lik,"Kaapsaoui") ~ "KAA",
                                                       str_detect(best_lik,"La_Grande") ~ "LAG",
                                                       str_detect(best_lik,"Beaver") ~ "BEA",
                                                       str_detect(best_lik,"Sabacunica") ~ "SAB",
                                                       str_detect(best_lik,"Eastmain") ~ "EAS",
                                                       str_detect(best_lik,"Rupert") ~ "RUP",
                                                       str_detect(best_lik,"Nottaway") ~ "NOT"))

whitefish_assignment %<>% mutate(site2 = case_when(str_detect(id,"ROG") ~ "ROG",
                                                  str_detect(id,"CHIm") ~ "KAA",
                                                  str_detect(id,"BIA") ~ "BIA",
                                                  str_detect(id,"LAG") & best_lik == "La_Grande" ~ "LAG",
                                                  str_detect(id,"LAG") & best_lik == "Roggan" ~ "ROG",
                                                  str_detect(id,"LAG") & best_lik == "Sabacunica" ~ "LAG",
                                                  str_detect(id,"BEA") ~ "BEA",
                                                  str_detect(id,"SAB") ~ "SAB",
                                                  str_detect(id,"MOA") ~ "SAB",
                                                  str_detect(id,"REN") ~ "REN",
                                                  str_detect(id,"TIL") ~ "REN",
                                                  str_detect(id,"CON") ~ "CON",
                                                  str_detect(id,"EAS") ~ "EAS",
                                                  str_detect(id,"FIS") ~ "EAS",
                                                  str_detect(id,"JAC") ~ "JAC",
                                                  str_detect(id,"RUP") ~ "RUP",
                                                  str_detect(id,"NOT") ~ "NOT"))



#Create a column which identifies individual as migrant

whitefish_assignment %<>% mutate(migrant = if_else(best_lik2 == site,"no","yes"))

whitefish_assignment %>% group_by(migrant) %>% dplyr::count(post_prob >= 0.8)

#--->
#Dispersal analysis

Nmax <- whitefish_assignment %>% filter(migrant=="yes",post_prob >= 0.8) %>% dplyr::count() %>% pull()

#Create vector of probabilities for random draw
probf <- as.vector(table(whitefish_assignment$site2)/sum(table(whitefish_assignment$site2)))

#expected random distribution of geographical distance 
#between source and misassigned rivers was generated 
#by drawing without replacement 10 000 pairs of populations

#rdraw <- lapply(1:10000, function(i) lapply(1:Nmax, function(i) sample(1:12, 2, replace = F, prob = probf)))

#saveRDS(rdraw, file = "data/dispersal_model/simulation_whitefish.rds")

rdraw <- read_rds("data/dispersal_model/simulation_whitefish.rds")

# Convert the list to a data frame
rdraw_df <- do.call(rbind, rdraw)

# Convert to a data frame and name the columns
rdraw_df <- as.data.frame(matrix(unlist(rdraw), ncol = 2, byrow = TRUE))
colnames(rdraw_df) <- c("pop1", "pop2")

# Add a column specifying the list number
rdraw_df$list_number <- rep(1:10000, each = Nmax)

# Add distance value

rdraw_dist <- left_join(rdraw_df, dfdist_wf, by = c("pop1","pop2"))

hist(rdraw_dist$values)

# Add distance class

rdraw_dist %<>% mutate(dist_class = case_when(values < 50 ~ "]0-50[",
                                              values >= 50 & values < 100 ~ "[50-100[",
                                              values >= 100 & values < 150 ~ "[100-150[",
                                              values >= 150 & values < 200 ~ "[150-200[",
                                              values >= 200 & values < 250 ~ "[200-250[",
                                              values >= 250 & values < 300 ~ "[250-300[",
                                              values >= 300 & values <= 355 ~ "[300-355["))

# Count number of individual falling in each class per list

rdraw_dist_count <- rdraw_dist %>% group_by(list_number) %>% dplyr::count(dist_class)

rdraw_dist_count$n <- rdraw_dist_count$n/Nmax

rdraw_dist_count %<>% group_by(dist_class) %>% dplyr::summarise(avg_n = mean(n),
                                                                lower_95ci = quantile(n, probs = 0.025),
                                                                upper_95ci = quantile(n,probs = 0.975))


#Relevel factors
rdraw_dist_count$dist_class <- as.factor(rdraw_dist_count$dist_class)

rdraw_dist_count$dist_class <- factor(rdraw_dist_count$dist_class,
                                      levels = c("]0-50[", "[50-100[", "[100-150[", "[150-200[", "[200-250[","[250-300[","[300-355["))


#Plot null distribution
ggplot(rdraw_dist_count) +
  geom_bar( aes(x=dist_class, y=avg_n), stat="identity", fill="white",color = "black", alpha=0.5) +
  geom_pointrange( aes(x=dist_class, y=avg_n, ymin=lower_95ci, ymax=upper_95ci), colour="darkred", alpha=0.9, size=1.3)+
  theme_bw()

#Average number of dispersers under random dispersal simulation

avg_dist <- rdraw_dist %>% group_by(list_number) %>% dplyr::summarise(mean_dist = mean(values))

hist(avg_dist$mean_dist)

mean(avg_dist$mean_dist)
quantile(avg_dist$mean_dist, probs = 0.025)
quantile(avg_dist$mean_dist, probs = 0.975)


# The average assignment distance of 81 individual whitefish randomly 
# dispersing over the 12 sampled populations by simulations was 
# 151 km [95%CI = (132 – 169)].

#--->
# Add observed distribution

site_name <- unique(whitefish_assignment$site2)
site_num <- as.numeric(as.factor(site_name))

df_name_num <- data.frame(site_name = site_name,
                          site_num = site_num)


whitefish_assignment %<>% left_join(df_name_num, by = c("site2"="site_name")) %>% dplyr::rename(pop1 = site_num)
whitefish_assignment %<>% left_join(df_name_num, by = c("best_lik2"="site_name")) %>% dplyr::rename(pop2 = site_num)

whitefish_disp <- whitefish_assignment %>% 
  filter(migrant=="yes",post_prob >= 0.8) %>% 
  dplyr::select(id,pop1,pop2) 

whitefish_disp %<>% left_join(dfdist_wf)

hist(log(whitefish_disp$values))

mean(whitefish_disp$values)
quantile(whitefish_disp$values, probs = 0.025)
quantile(whitefish_disp$values, probs = 0.975)

t.test(log(whitefish_disp$values), log(avg_dist$mean_dist))

#In whitefish, the 81 misassigned individuals were assigned at 
#an average distance of CI(12 - 345 km from their sampling river

whitefish_disp  %<>% mutate(dist_class = case_when(values < 50 ~ "]0-50[",
                                               values >= 50 & values < 100 ~ "[50-100[",
                                               values >= 100 & values < 150 ~ "[100-150[",
                                               values >= 150 & values < 200 ~ "[150-200[",
                                               values >= 200 & values < 250 ~ "[200-250[",
                                               values >= 250 & values < 300 ~ "[250-300[",
                                               values >= 300 & values <= 355 ~ "[300-355["))

df_disp_obs <- whitefish_disp %>% dplyr::count(dist_class)

df_disp_obs$n <- df_disp_obs$n/Nmax

df_disp_obs <- rbind(df_disp_obs,c("[200-250[",0))
df_disp_obs <- rbind(df_disp_obs,c("[250-300[",0))

df_disp_obs$dist_class <- as.factor(df_disp_obs$dist_class)

df_disp_obs$n <- as.numeric(df_disp_obs$n)

df_disp_obs$dist_class <- factor(df_disp_obs$dist_class,
                                 levels = c("]0-50[", "[50-100[", "[100-150[", "[150-200[", "[200-250[","[250-300[","[300-355["))

ggplot() +
  geom_ribbon(data = rdraw_dist_count,aes(ymin = lower_95ci, ymax = upper_95ci,x= dist_class,group=1), fill = "grey70",alpha = 0.4,colour = "black",linetype=2) +
  geom_line(data = rdraw_dist_count,aes(y = avg_n, x = dist_class,group=1))+
  geom_point(data = df_disp_obs, aes(x = dist_class, y = n),color="darkred")+
  geom_line(data = df_disp_obs, aes(x = dist_class, y = n,group=1),color="darkred")+
  ylab("Proportion of dispersers")+
  xlab("\nDispersal distance class (km)")+
  ylim(c(0,1))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle = 45,vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("figures/dispersal_whitefish.png",width =7, height = 7)


rdraw_dist_count


#GTEST whitefish

gt_disp_obs <- whitefish_disp %>% dplyr::count(dist_class)

gt_disp_obs <- rbind(gt_disp_obs,c("[200-250[",0))
gt_disp_obs <- rbind(gt_disp_obs,c("[250-300[",0))

rdraw_dist_count$n <- rdraw_dist_count$avg_n*81

gt_wf <- gt_disp_obs %>% left_join(rdraw_dist_count, by = "dist_class")

gt_wf_f <- rbind(as.numeric(gt_wf$n.x),gt_wf$n.y)

GTest(gt_wf_f, correct = "none")

#Gt test between species


sp_gt <- rbind(gt_bt_f[1,],gt_wf_f[1,])

GTest(sp_gt, correct = "none")

