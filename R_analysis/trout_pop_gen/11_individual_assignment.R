rm(list = ls())

#-----------------#
#---- LIBRARY ----#
#-----------------#

library(tidyverse)
library(marmap)

#-----------------#
#----- DATA ------#
#-----------------#

df <- read.table("data/individual_assignment_analysis/ind_assignment.txt")

df$pop1 <- as.numeric(as.factor(df$pop))

df %<>% mutate(pop2_tmp = case_when(cluster == "KAA" & pop %in% c("KAA","ROG") ~ pop,
                                    T ~ cluster))

df$pop2 <- as.numeric(as.factor(df$pop2_tmp))

df <- left_join(df, dfdist, by = c("pop1","pop2"))

nrow(df[df$values > 0,])

#There is a total of 158 dispersers in the dataset

#-----------------#
#---- ANALYSIS ---#
#-----------------#

#--->

JB <- getNOAA.bathy(lon1 = -82.084, lon2 = -78,
                    lat1 = 51.021, lat2 = 55, resolution = 1)

sites <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/pop_sampled.txt",header=T)

sites_bt <- plyr::arrange(sites, pop_code, desc(pop_code)) %>% filter(!pop_code %in% c("RUP","NOT","JAC"))

marsites <- sites_bt %>% dplyr::select(long,lat) %>% dplyr::rename(x=long,y=lat)

trans1 <- trans.mat(JB,min.depth = 5)

dist1 <- lc.dist(trans1,marsites,res="dist")

dist1

mat_dist<-as.matrix(dist1)

dfdist <- setNames(melt(mat_dist), c('pop1', 'pop2', 'values'))

# Join with dataset to get number of dispersers

df <- left_join(df, dfdist, by = c("pop1","pop2"))

nrow(df[df$values > 0,])

#There is a total of 75 dispersers in the dataset

#--->

Nmax <- nrow(df[df$values > 0,])

#Create vector of probabilities for random draw
probf <- as.vector(table(df$pop)/sum(table(df$pop)))

#expected random distribution of geographical distance 
#between source and misassigned rivers was generated 
#by drawing without replacement 10 000 pairs of populations

#rdraw <- lapply(1:10000, function(i) lapply(1:Nmax, function(i) sample(1:11, 2, replace = F, prob = probf)))

saveRDS(rdraw, file = "data/individual_assignment_analysis/simulation.rds")

rdraw <- read_rds("data/individual_assignment_analysis/simulation.rds")

# Convert the list to a data frame
rdraw_df <- do.call(rbind, rdraw)

# Convert to a data frame and name the columns
rdraw_df <- as.data.frame(matrix(unlist(rdraw), ncol = 2, byrow = TRUE))
colnames(rdraw_df) <- c("pop1", "pop2")

# Add a column specifying the list number
rdraw_df$list_number <- rep(1:10000, each = Nmax)

# Add distance value

rdraw_dist <- left_join(rdraw_df, dfdist, by = c("pop1","pop2"))

#hist(rdraw_dist$values)

# Add distance class

rdraw_dist %<>% mutate(dist_class = case_when(values < 50 ~ "]0-50[",
                                              values >= 50 & values < 100 ~ "[50-100[",
                                              values >= 100 & values < 150 ~ "[100-150[",
                                              values >= 150 & values < 200 ~ "[150-200[",
                                              values >= 200 & values < 250 ~ "[200-250[",
                                              values >= 250 & values < 300 ~ "[250-300[",
                                              values >= 300 & values <= 355 ~ "[300 - 355["))

# Count number of individual falling in each class per list

rdraw_dist_count <- rdraw_dist %>% group_by(list_number) %>% dplyr::count(dist_class)

rdraw_dist_count$n <- rdraw_dist_count$n/Nmax

rdraw_dist_count %<>% group_by(dist_class) %>% dplyr::summarise(avg_n = mean(n),
                                                                lower_95ci = quantile(n, probs = 0.025),
                                                                upper_95ci = quantile(n,probs = 0.975))


#Relevel factors
rdraw_dist_count$dist_class <- as.factor(rdraw_dist_count$dist_class)

rdraw_dist_count$dist_class <- factor(rdraw_dist_count$dist_class,
                                      levels = c("]0-50[", "[50-100[", "[100-150[", "[150-200[", "[200-250[","[250-300[","[300 - 355["))


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

# The average assignment distance of 75 individual brook trout randomly 
# dispersing over the 11 sampled populations by simulations was 
# 131 km [95%CI = (113.5 – 149.9)].


#--->
# Add observed distribution

df_disp <- df %>% filter(values > 0)

hist(df_disp$values)

mean(df_disp$values)

#In brook charr, the 75 misassigned individuals were assigned at 
#an average distance of 49 km from their sampling river

df_disp %<>% mutate(dist_class = case_when(values < 50 ~ "]0-50[",
                                           values >= 50 & values < 100 ~ "[50-100[",
                                           values >= 100 & values < 150 ~ "[100-150[",
                                           values >= 150 & values < 200 ~ "[150-200[",
                                           values >= 200 & values < 250 ~ "[200-250[",
                                           values >= 250 & values < 300 ~ "[250-300[",
                                           values >= 300 & values <= 355 ~ "[300 - 355["))

df_disp_obs <- df_disp %>% dplyr::count(dist_class)

df_disp_obs$n <- df_disp_obs$n/Nmax

df_disp_obs$dist_class <- as.factor(df_disp_obs$dist_class)

df_disp_obs$dist_class <- factor(df_disp_obs$dist_class,
                                 levels = c("]0-50[", "[50-100[", "[100-150[", "[150-200[", "[200-250[","[250-300[","[300 - 355["))

to_add <- data.frame(dist_class = c("[150-200[","[200-250[","[250-300[","[300 - 355["),
           n = c(0,0,0,0))

df_disp_obs <- rbind(df_disp_obs,to_add)

df_disp_obs$n <- as.numeric(df_disp_obs$n)

#Plot null distribution + observed data
ggplot() +
  geom_bar(data = rdraw_dist_count,aes(x=dist_class, y=avg_n), stat="identity", fill="white",color = "black", alpha=0.5) +
  geom_point(data = df_disp_obs, aes(x = dist_class, y = n))+
  geom_line(data = df_disp_obs, aes(x = dist_class, y = n,group=1))+
  geom_pointrange(data = rdraw_dist_count, aes(x=dist_class, y=avg_n, ymin=lower_95ci, ymax=upper_95ci), colour="darkred", alpha=0.9, size=1.3)+
  ylab("Proportion of dispersers")+
  xlab("Dispersal distance class (km)")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle = 45,vjust = 0.5),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("figures/obs_disp_distance_vs_exp.jpeg")

saveRDS(rdraw_dist_count,"~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/individual_assignment_analysis/random_disp_sim_bt.rds")

saveRDS(df_disp_obs,"~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/individual_assignment_analysis/obs_disp_bt.rds")


# G-test to compare distribution

obs_exp <- left_join(df_disp_obs, rdraw_dist_count, by = "dist_class")

saveRDS(obs_exp,"~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/individual_assignment_analysis/trout_misassigned.rds")

obs_exp_mat <- rbind(obs_exp$n * Nmax,obs_exp$avg_n * Nmax)

DescTools::GTest(obs_exp_mat,correct = "none")






