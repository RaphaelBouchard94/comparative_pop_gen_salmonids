rm(list = ls())

####################################
###########Library##################
####################################

library(tidyverse)
library(magrittr)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(patchwork)

####################################
######  Admixture analysis  ########
####################################

#import ngsadmix proportion of admixture

files <- list.files(path = "~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/ngs_admix", 
                    pattern = "qopt",
                    full.names = T)


for (i in 1:length(files)){
  new_df <-  read.table(files[i])
  assign(substr(files[i],101,108),new_df)
}

#load samples_id in short format
id <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/ind_list_without_missing.txt")


id[duplicated(id$V1),]

colnames(id)<-"id"

#check if duplicated if
id[duplicated(id$V1),]

#row_numer col
id %<>% dplyr::mutate(rown = as.factor(row_number()),
               pop = case_when(str_detect(string = id, "BEA") ~ "BEA",
                               str_detect(string = id, "FIS") ~ "EAS",
                               str_detect(string = id, "JAC") ~ "JAC",
                               str_detect(string = id, "SAB") ~ "SAB",
                               str_detect(string = id, "BIA") ~ "BIA",
                               str_detect(string = id, "LAG") ~ "LAG",
                               str_detect(string = id, "TIL") ~ "REN",
                               str_detect(string = id, "MOA") ~ "SAB",
                               str_detect(string = id, "ROG") ~ "ROG",
                               str_detect(string = id, "EAS") ~ "EAS",
                               str_detect(string = id, "CHIm") ~ "KAA",
                               str_detect(string = id, "CON") ~ "CON",
                               str_detect(string = id, "AQU") ~ "AQU",
                               str_detect(string = id, "REN") ~ "REN",
                               str_detect(string = id, "BRO") ~ "BRO",
                               str_detect(string = id, "WAS") ~ "WAS"))
                                     

id$pop<-factor(id$pop, 
               levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","EAS","BRO"))

####################################
#Check best cluster

log_lik <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/ngs_admix/k_lik.txt")

colnames(log_lik) <- "likelihood"

rep_lik <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/ngs_admix/k_rep.txt")

colnames(rep_lik) <- c("K","rep")

lik_df <- as.data.frame(cbind(log_lik$likelihood,rep_lik$K,rep_lik$rep))

colnames(lik_df) <- c("likelihood","K","rep")

lik_df %>%  ggplot(aes(x=K,y=likelihood))+
  geom_point()+
  theme_bw()

logk <- tapply(lik_df$likelihood, lik_df$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

bestk$K_reorder <- fct_relevel(bestk$K,c("K1","K2","K3","K4","K5","K6","K7","K8","K9","K10","K11"))

bestk <- data.frame(loglik = logk,K = rownames(logk))

bestk %>% 
  filter(!K == 1) %>% 
  ggplot(aes(x=as.numeric(K),y=loglik))+
  geom_point()+
  geom_line()+
  scale_x_continuous(breaks=seq(1, 14, 1))+
  xlab("")+
  ylab("logLik")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/loglik_K.jpeg",height=5,width = 6)


lik_df %>% filter(K != 1) %>% group_by(as.factor(K)) %>% top_n(1, likelihood)

####################################

#---> Start by making dataframes for the relevant structure plot 

trout_apclus <- read_rds("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/ap_cluster/ap_cluster_trout.Rds")

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

####################################
#Palette


pal_raph <- c("#abdda4", #Beaver
                 "#F67E4BFF", #Eastmain
                 "#4575b4", #LaGrande
                 "#d53e4f", #Broadback
                 "#74add1", #Bia
                 "#FEDA8BFF", #Renoyer
                 "#66c2a5", #Aquatuc
                 "#FDB366FF", #Conn
                 "#313695", #Roggan/Kaapsaoui
                 "#EAECCCFF",
              "#607D8BFF",
              "#8BC34AFF",
              "#795548FF")


####################################
#K2

k2 <- K2.qopt %>% 
  dplyr::mutate(id = as.integer(row_number())) %>% 
  gather('pop', 'prob', V1:V2) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k2 %<>% left_join(id,by=c("id"="rown"))

k2 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k2 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                      str_detect(id.y,"CHIm") ~ "KAA",
                                      str_detect(id.y,"BIA") ~ "BIA",
                                      str_detect(id.y,"LAG") ~ "LAG",
                                      str_detect(id.y,"AQU") ~ "AQU",
                                      str_detect(id.y,"BEA") ~ "BEA",
                                      str_detect(id.y,"SAB") ~ "SAB",
                                      str_detect(id.y,"MOA") ~ "SAB",
                                      str_detect(id.y,"REN") ~ "REN",
                                      str_detect(id.y,"CON") ~ "CON",
                                      str_detect(id.y,"EAS") ~ "EAS",
                                      str_detect(id.y,"FIS") ~ "EAS",
                                      str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                      str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                      str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                      str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                      str_detect(id.y,"RUP") ~ "RUP",
                                      str_detect(id.y,"BRO") ~ "BRO",
                                      str_detect(id.y,"NOT") ~ "NOT"))

k2$site<-factor(k2$site, 
                 levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


a<-ggplot(k2, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K2")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=18),
        legend.position = "none")

####################################
#K3

k3 <-K3.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V3) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k3 %<>% left_join(id,by=c("id"="rown"))

k3 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k3 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                str_detect(id.y,"CHIm") ~ "KAA",
                                str_detect(id.y,"BIA") ~ "BIA",
                                str_detect(id.y,"LAG") ~ "LAG",
                                str_detect(id.y,"AQU") ~ "AQU",
                                str_detect(id.y,"BEA") ~ "BEA",
                                str_detect(id.y,"SAB") ~ "SAB",
                                str_detect(id.y,"MOA") ~ "SAB",
                                str_detect(id.y,"REN") ~ "REN",
                                str_detect(id.y,"CON") ~ "CON",
                                str_detect(id.y,"EAS") ~ "EAS",
                                str_detect(id.y,"FIS") ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                str_detect(id.y,"RUP") ~ "RUP",
                                str_detect(id.y,"BRO") ~ "BRO",
                                str_detect(id.y,"NOT") ~ "NOT"))

k3$site<-factor(k3$site, 
                levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


b<-ggplot(k3, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K3")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")



####################################
#K4

k4 <- K4.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k4 %<>% left_join(id,by=c("id"="rown"))

k4 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k4 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                str_detect(id.y,"CHIm") ~ "KAA",
                                str_detect(id.y,"BIA") ~ "BIA",
                                str_detect(id.y,"LAG") ~ "LAG",
                                str_detect(id.y,"AQU") ~ "AQU",
                                str_detect(id.y,"BEA") ~ "BEA",
                                str_detect(id.y,"SAB") ~ "SAB",
                                str_detect(id.y,"MOA") ~ "SAB",
                                str_detect(id.y,"REN") ~ "REN",
                                str_detect(id.y,"CON") ~ "CON",
                                str_detect(id.y,"EAS") ~ "EAS",
                                str_detect(id.y,"FIS") ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                str_detect(id.y,"RUP") ~ "RUP",
                                str_detect(id.y,"BRO") ~ "BRO",
                                str_detect(id.y,"NOT") ~ "NOT"))

k4$site<-factor(k4$site, 
                levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


c<-ggplot(k4, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K4")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

####################################
#K5

k5 <- K5.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V5) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k5 %<>% left_join(id,by=c("id"="rown"))

k5 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k5 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                str_detect(id.y,"CHIm") ~ "KAA",
                                str_detect(id.y,"BIA") ~ "BIA",
                                str_detect(id.y,"LAG") ~ "LAG",
                                str_detect(id.y,"AQU") ~ "AQU",
                                str_detect(id.y,"BEA") ~ "BEA",
                                str_detect(id.y,"SAB") ~ "SAB",
                                str_detect(id.y,"MOA") ~ "SAB",
                                str_detect(id.y,"REN") ~ "REN",
                                str_detect(id.y,"CON") ~ "CON",
                                str_detect(id.y,"EAS") ~ "EAS",
                                str_detect(id.y,"FIS") ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                str_detect(id.y,"RUP") ~ "RUP",
                                str_detect(id.y,"BRO") ~ "BRO",
                                str_detect(id.y,"NOT") ~ "NOT"))

k5$site<-factor(k5$site, 
                levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


d<-ggplot(k5, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K5")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

####################################
#K6

k6 <- K6.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V6) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k6 %<>% left_join(id,by=c("id"="rown"))

k6 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k6 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                str_detect(id.y,"CHIm") ~ "KAA",
                                str_detect(id.y,"BIA") ~ "BIA",
                                str_detect(id.y,"LAG") ~ "LAG",
                                str_detect(id.y,"AQU") ~ "AQU",
                                str_detect(id.y,"BEA") ~ "BEA",
                                str_detect(id.y,"SAB") ~ "SAB",
                                str_detect(id.y,"MOA") ~ "SAB",
                                str_detect(id.y,"REN") ~ "REN",
                                str_detect(id.y,"CON") ~ "CON",
                                str_detect(id.y,"EAS") ~ "EAS",
                                str_detect(id.y,"FIS") ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                str_detect(id.y,"RUP") ~ "RUP",
                                str_detect(id.y,"BRO") ~ "BRO",
                                str_detect(id.y,"NOT") ~ "NOT"))

k6$site<-factor(k6$site, 
                levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


e<-ggplot(k6, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K6")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")


####################################
#K7

k7 <- K7.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V7) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k7 %<>% left_join(id,by=c("id"="rown"))

k7 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k7 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                str_detect(id.y,"CHIm") ~ "KAA",
                                str_detect(id.y,"BIA") ~ "BIA",
                                str_detect(id.y,"LAG") ~ "LAG",
                                str_detect(id.y,"AQU") ~ "AQU",
                                str_detect(id.y,"BEA") ~ "BEA",
                                str_detect(id.y,"SAB") ~ "SAB",
                                str_detect(id.y,"MOA") ~ "SAB",
                                str_detect(id.y,"REN") ~ "REN",
                                str_detect(id.y,"CON") ~ "CON",
                                str_detect(id.y,"EAS") ~ "EAS",
                                str_detect(id.y,"FIS") ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                str_detect(id.y,"RUP") ~ "RUP",
                                str_detect(id.y,"BRO") ~ "BRO",
                                str_detect(id.y,"NOT") ~ "NOT"))

k7$site<-factor(k7$site, 
                levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))



f<-ggplot(k7, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K7")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

####################################
#K8

k8 <- K8.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V8) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k8 %<>% left_join(id,by=c("id"="rown"))

k8 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k8 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                str_detect(id.y,"CHIm") ~ "KAA",
                                str_detect(id.y,"BIA") ~ "BIA",
                                str_detect(id.y,"LAG") ~ "LAG",
                                str_detect(id.y,"AQU") ~ "AQU",
                                str_detect(id.y,"BEA") ~ "BEA",
                                str_detect(id.y,"SAB") ~ "SAB",
                                str_detect(id.y,"MOA") ~ "SAB",
                                str_detect(id.y,"REN") ~ "REN",
                                str_detect(id.y,"CON") ~ "CON",
                                str_detect(id.y,"EAS") ~ "EAS",
                                str_detect(id.y,"FIS") ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                str_detect(id.y,"RUP") ~ "RUP",
                                str_detect(id.y,"BRO") ~ "BRO",
                                str_detect(id.y,"NOT") ~ "NOT"))

k8$site<-factor(k8$site, 
                levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))



g<-ggplot(k8, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K8")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

####################################
#K9

k9 <- K9.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V9) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k9 %<>% left_join(id,by=c("id"="rown"))
                  
k9 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k9 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                str_detect(id.y,"CHIm") ~ "KAA",
                                str_detect(id.y,"BIA") ~ "BIA",
                                str_detect(id.y,"LAG") ~ "LAG",
                                str_detect(id.y,"AQU") ~ "AQU",
                                str_detect(id.y,"BEA") ~ "BEA",
                                str_detect(id.y,"SAB") ~ "SAB",
                                str_detect(id.y,"MOA") ~ "SAB",
                                str_detect(id.y,"REN") ~ "REN",
                                str_detect(id.y,"CON") ~ "CON",
                                str_detect(id.y,"EAS") ~ "EAS",
                                str_detect(id.y,"FIS") ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                str_detect(id.y,"RUP") ~ "RUP",
                                str_detect(id.y,"BRO") ~ "BRO",
                                str_detect(id.y,"NOT") ~ "NOT"))

k9$site<-factor(k9$site, 
                levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


h<-ggplot(k9, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K9")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")


####################################
#K10

k10 <- K10.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V10) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k10 %<>% left_join(id,by=c("id"="rown"))

k10 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k10 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                str_detect(id.y,"CHIm") ~ "KAA",
                                str_detect(id.y,"BIA") ~ "BIA",
                                str_detect(id.y,"LAG") ~ "LAG",
                                str_detect(id.y,"AQU") ~ "AQU",
                                str_detect(id.y,"BEA") ~ "BEA",
                                str_detect(id.y,"SAB") ~ "SAB",
                                str_detect(id.y,"MOA") ~ "SAB",
                                str_detect(id.y,"REN") ~ "REN",
                                str_detect(id.y,"CON") ~ "CON",
                                str_detect(id.y,"EAS") ~ "EAS",
                                str_detect(id.y,"FIS") ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                str_detect(id.y,"RUP") ~ "RUP",
                                str_detect(id.y,"BRO") ~ "BRO",
                                str_detect(id.y,"NOT") ~ "NOT"))

k10$site<-factor(k10$site, 
                levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


i<-ggplot(k10, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K10")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")


####################################
#K11

k11 <- K11.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V11) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k11 %<>% left_join(id,by=c("id"="rown"))

k11 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k11 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                 str_detect(id.y,"CHIm") ~ "KAA",
                                 str_detect(id.y,"BIA") ~ "BIA",
                                 str_detect(id.y,"LAG") ~ "LAG",
                                 str_detect(id.y,"AQU") ~ "AQU",
                                 str_detect(id.y,"BEA") ~ "BEA",
                                 str_detect(id.y,"SAB") ~ "SAB",
                                 str_detect(id.y,"MOA") ~ "SAB",
                                 str_detect(id.y,"REN") ~ "REN",
                                 str_detect(id.y,"CON") ~ "CON",
                                 str_detect(id.y,"EAS") ~ "EAS",
                                 str_detect(id.y,"FIS") ~ "EAS",
                                 str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                 str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                 str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                 str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                 str_detect(id.y,"RUP") ~ "RUP",
                                 str_detect(id.y,"BRO") ~ "BRO",
                                 str_detect(id.y,"NOT") ~ "NOT"))

k11$site<-factor(k11$site, 
                 levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


j<-ggplot(k11, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K11")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

####################################
#K12

k12 <- K12.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V12) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k12 %<>% left_join(id,by=c("id"="rown"))

k12 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k12 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                 str_detect(id.y,"CHIm") ~ "KAA",
                                 str_detect(id.y,"BIA") ~ "BIA",
                                 str_detect(id.y,"LAG") ~ "LAG",
                                 str_detect(id.y,"AQU") ~ "AQU",
                                 str_detect(id.y,"BEA") ~ "BEA",
                                 str_detect(id.y,"SAB") ~ "SAB",
                                 str_detect(id.y,"MOA") ~ "SAB",
                                 str_detect(id.y,"REN") ~ "REN",
                                 str_detect(id.y,"CON") ~ "CON",
                                 str_detect(id.y,"EAS") ~ "EAS",
                                 str_detect(id.y,"FIS") ~ "EAS",
                                 str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                 str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                 str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                 str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                 str_detect(id.y,"RUP") ~ "RUP",
                                 str_detect(id.y,"BRO") ~ "BRO",
                                 str_detect(id.y,"NOT") ~ "NOT"))

k12$site<-factor(k12$site, 
                 levels = c("ROG","KAA","BIA","LAG","AQU","BEA","MOA","SAB","TIL","REN","CON","FIS","EAS","JAC","RUP","BRO","NOT"))


k<-ggplot(k12, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K12")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")


####################################
#K13 

k13 <- K13.qopt %>% 
  dplyr::mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V13) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k13 %<>% left_join(id,by=c("id"="rown"))

k13 %<>% left_join(trout_apclus,by=c("id.y"="id"))

k13 %<>% mutate(site = case_when(str_detect(id.y,"ROG") ~ "ROG",
                                 str_detect(id.y,"CHIm") ~ "KAA",
                                 str_detect(id.y,"BIA") ~ "BIA",
                                 str_detect(id.y,"LAG") ~ "LAG",
                                 str_detect(id.y,"AQU") ~ "AQU",
                                 str_detect(id.y,"BEA") ~ "BEA",
                                 str_detect(id.y,"SAB") ~ "SAB",
                                 str_detect(id.y,"MOA") ~ "SAB",
                                 str_detect(id.y,"REN") ~ "REN",
                                 str_detect(id.y,"CON") ~ "CON",
                                 str_detect(id.y,"EAS") ~ "EAS",
                                 str_detect(id.y,"FIS") ~ "EAS",
                                 str_detect(id.y,"WASm") & cluster == "BRO" ~ "BRO",
                                 str_detect(id.y,"WASm") & cluster == "EAS" ~ "EAS",
                                 str_detect(id.y,"WASm") & cluster == "CON" ~ "CON",
                                 str_detect(id.y,"WASm") & cluster == "REN" ~ "REN",
                                 str_detect(id.y,"RUP") ~ "RUP",
                                 str_detect(id.y,"BRO") ~ "BRO",
                                 str_detect(id.y,"NOT") ~ "NOT"))


l<-ggplot(k13, aes(id, prob, fill = pop.x)) +
  geom_col(width=2) +
  scale_fill_manual(values = pal_raph)+
  facet_grid(~site,scales = "free_x")+
  ylab("K13")+
  xlab("")+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")


####################################

a/b/c/d/e/f/g/h/i/j/k/l

ggsave("figures/k2_k13.jpeg",height=15,width = 20)



