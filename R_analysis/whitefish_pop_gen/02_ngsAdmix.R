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

files <- list.files(path = "~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/ngs_admix", 
                    pattern = "qopt",
                    full.names = T)

for (i in 1:length(files)){
  new_df <-  read.table(files[i])
  assign(paste("k",substr(files[i],156,159),sep=""),new_df)
}

#load samples_id in short format
id <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/ngs_admix/id_ind_bam.tsv")

colnames(id)<-"id"

#check if duplicated if
id[duplicated(id$V1),]

#row_numer col
id %<>% mutate(rown = as.factor(row_number()),
               pop = case_when(str_detect(string = id, "RUP") ~ "RUP",
                               str_detect(string = id, "NOT") ~ "NOT",
                               str_detect(string = id, "BEA") ~ "BEA",
                               str_detect(string = id, "FIS") ~ "FIS",
                               str_detect(string = id, "JAC") ~ "JAC",
                               str_detect(string = id, "SAB") ~ "SAB",
                               str_detect(string = id, "BIA") ~ "BIA",
                               str_detect(string = id, "LAG") ~ "LAG",
                               str_detect(string = id, "TIL") ~ "TIL",
                               str_detect(string = id, "MOA") ~ "MOA",
                               str_detect(string = id, "ROG") ~ "ROG",
                               str_detect(string = id, "EAS") ~ "EAS",
                               str_detect(string = id, "CHIm") ~ "KAA",
                               str_detect(string = id, "CON") ~ "CON"))

id$pop<-factor(id$pop, 
                  levels = c("ROG","KAA","BIA","LAG","BEA","SAB","TIL","CON","FIS","EAS","JAC","RUP","NOT"))

####################################
#Check best cluster

log_lik <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/ngs_admix/k_lik.txt")

colnames(log_lik) <- "likelihood"

rep_lik <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/ngs_admix/k_rep.txt")

colnames(rep_lik) <- c("K","rep")

lik_df <- as.data.frame(cbind(log_lik$likelihood,rep_lik$K,rep_lik$rep))

colnames(lik_df) <- c("likelihood","K","rep")

lik_df %>%  ggplot(aes(x=K,y=likelihood))+
  geom_point()+
  theme_bw()

logk <- tapply(lik_df$likelihood, lik_df$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

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

ggsave("figure/loglik_K.jpeg",height=5,width = 6)

ggsave("~/Desktop/figures_ROC_june/loglik_K.jpeg",height=5,width = 6)

####################################
#Palette

pal_raph<-c("#2196F3FF",
            "#C62828FF",
            "#FF9800FF",
            "#4CAF50FF",
            "#006064FF",
            "#9C27B0FF",
            "#607D8BFF",
            "#8BC34AFF",
            "#795548FF",
            "#3F51B5FF",
            "#607D8BFF",
            "#293352", 
            "#C4961A")


####################################
#K2

k2 <- k_2_1 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V2) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k2 %<>% left_join(id,by=c("id"="rown"))

a<-ggplot(k2, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k3 <- k_3_1 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V3) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k3 %<>% left_join(id,by=c("id"="rown"))

b<-ggplot(k3, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k4 <- k_4_1 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V4) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k4 %<>% left_join(id,by=c("id"="rown"))

c<-ggplot(k4, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k5 <- k_5_1 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V5) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k5 %<>% left_join(id,by=c("id"="rown"))

d<-ggplot(k5, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k6 <- k_6_1 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V6) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k6 %<>% left_join(id,by=c("id"="rown"))

e<-ggplot(k6, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k7 <- k_7_1 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V7) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k7 %<>% left_join(id,by=c("id"="rown"))

f<-ggplot(k7, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k8 <- k_8_1 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V8) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k8 %<>% left_join(id,by=c("id"="rown"))

g<-ggplot(k8, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k9 <- k_9_1 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V9) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k9 %<>% left_join(id,by=c("id"="rown"))

h<-ggplot(k9, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k10 <- k_10_ %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V10) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k10 %<>% left_join(id,by=c("id"="rown"))

i<-ggplot(k10, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k11 <- k_11_ %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V11) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k11 %<>% left_join(id,by=c("id"="rown"))

j<-ggplot(k11, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k12 <- k_12_ %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V12) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k12 %<>% left_join(id,by=c("id"="rown"))

k<-ggplot(k12, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

k13 <- k_13_ %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V13) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k13 %<>% left_join(id,by=c("id"="rown"))

l <- ggplot(k13, aes(id, prob, fill = pop.x)) +
  geom_col() +
  scale_fill_manual(values=pal_raph)+
  theme_classic()+
  facet_grid(~pop.y,scales = "free_x")+
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

########

patchwork<-a/b/c/d/e/f/g/h/i/j/k/l

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/figure/k2_k13.jpeg",height=15,width = 20)
  









