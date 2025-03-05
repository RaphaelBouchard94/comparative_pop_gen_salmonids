rm(list = ls())

#####################
#######Library#######
#####################

library(tidyverse)
library(magrittr)
library(paletteer)
library(plotly)

#####################
#########DATA########
#####################

#####################
#Based on PCA
cov <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/all_maf0.05pctind0.85_maxdepth10_ALL_CHR.pruned_maxdist_200kb_weight_0125_singletons.pcangsd.cov")

#Make PCA using covariance matrix
pca <- eigen(cov)
pca.mat<-as.matrix(pca$vectors %*% (diag(pca$values))^0.5)

#add column names
nPC<-dim(pca$vectors)[2]
col_PC<-vector(length=nPC)
for (i in 1 : nPC) {col_PC[i]<-paste0("PC",i)}
colnames(pca.mat)<-c(col_PC)

#add rownames
samples <- read.table("id_ind_bam.tsv")
rownames(pca.mat)<-samples$V1

#calculate varsum(eigen_mats$values[eigen_mats$values>=0]
var1<-round(pca$values[1]*100/sum(pca$values[pca$values>=0]),2)
var2<-round(pca$values[2]*100/sum(pca$values[pca$values>=0]),2)
var3<-round(pca$values[3]*100/sum(pca$values[pca$values>=0]),2)
var4<-round(pca$values[4]*100/sum(pca$values[pca$values>=0]),2)


#Make dataset to plot
pca_cocl <- data.frame(pca.mat[,1:4])

#sample id column
pca_cocl$id <- samples$V1

#pop column
pca_cocl %<>% mutate(pop = case_when(str_detect(string = id, "RUP") ~ "RUP",
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

##add info on pop
pop_file <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/pop_sampled.txt",header=T)

pop_file %>% arrange(desc(lat)) 

pca_cocl %<>% left_join(pop_file,by = c("pop" = "pop_code"))

pca_cocl %>% group_by(pop) %>% count()

ggplot(pca_cocl, aes(x = PC1,y = PC2,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC1 (3.81%)")+
  ylab("PC2 (1.19%)")+
  scale_color_manual(values = paletteer_d(`"awtools::bpalette"`), name = "Population (n)")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


kaa <-pca_cocl %>% filter(pop == "KAA",PC2>0.1) %>% select(id) %>% mutate(genetic_cluster = rep("KAA",26))

rog <- pca_cocl %>% filter(PC1>0.28, PC2 <0.1) %>% select(id) %>% mutate(genetic_cluster = rep("ROG",35))

lag_bia <- pca_cocl %>% filter(PC1>0.09, PC2 <0)%>% select(id) %>% mutate(genetic_cluster = rep("LAG",45))

eas <- pca_cocl %>% filter(PC1< -0.08,PC2 <0,PC2>-0.07) %>% select(id) %>% mutate(genetic_cluster = rep("EAS",77))

rup <-pca_cocl %>% filter(PC1< -0.125,PC2 < 1.72,PC2>0.04) %>% select(id) %>% mutate(genetic_cluster = rep("RUP",105))

not<- pca_cocl %>% filter(PC1< -0.125,PC2>0.2) %>% select(id) %>% mutate(genetic_cluster = rep("NOT",15))

#####################
#Based on K13
#load samples_id in short format
id <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/ngs_admix/id_ind_bam.tsv")

colnames(id)<-"id"

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


k13 <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/ngs_admix/all_maf0.05_pctind0.85_maxdepth10_pruned_singletons_13_1.qopt")

k13 <- k13 %>% 
  mutate(id = row_number()) %>% 
  gather('pop', 'prob', V1:V13) %>% 
  group_by(id) %>% 
  mutate(likely_assignment = pop[which.max(prob)],
         assingment_prob = max(prob)) %>% 
  arrange(likely_assignment, desc(assingment_prob)) %>% 
  ungroup() %>% 
  mutate(id = forcats::fct_inorder(factor(id)))

k13 %<>% left_join(id,by=c("id"="rown"))

ggplot(k13, aes(id, prob, fill = pop.x)) +
  geom_col() +
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
        axis.ticks.x = element_blank())


sab <- k13 %>% filter(pop.y %in% c("SAB","BEA"),likely_assignment == "V10"& prob > 0.5) %>% select(id.y) %>% rename(id = id.y) %>% mutate(genetic_cluster = rep("SAB",41))

bea <- k13 %>% filter(pop.y %in% c("SAB","BEA"),likely_assignment == "V9"& prob > 0.5)%>% select(id.y) %>% rename(id = id.y) %>% mutate(genetic_cluster = rep("BEA",55))


##############
#Final dataset

genetic_cluster<-bind_rows(rog,kaa,lag_bia,bea,sab,eas,rup,not)

write.table(genetic_cluster,"~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/distinct_genetic_cluster/true_genet_cluster.txt",
            quote = F,
            row.names = F,
            col.names = F)










