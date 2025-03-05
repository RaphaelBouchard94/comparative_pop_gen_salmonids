rm(list = ls())

#####################
#######Library#######
#####################

library(tidyverse)
library(magrittr)
library(paletteer)
library(plotly)
library(MASS)
library(ggridges)
library(caret)
library(klaR)

#####################
#########DATA########
#####################

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

#add true genetic cluster column

pop <- read.table("data/distinct_genetic_cluster/true_genet_cluster.txt")

dup <- pop[which(duplicated(pop$V1)),]

not <- pop %>% filter(V1 %in% dup$V1,V2=="NOT")

remove <- pop %>% filter(V1 %in% dup$V1)

pop %<>% filter(!V1 %in% remove$V1) %>% bind_rows(not)

pca_cocl %<>% left_join(pop,by=c("id"="V1")) %>% rename(pop=V2)

##add info on pop
pop_file <- read.table("data/fst/pop_sampled_whitefish.txt",header=T)

pop_file %>% arrange(desc(lat)) 

pca_cocl %<>% left_join(pop_file,by = c("pop" = "pop_code"))

pca_cocl %>% group_by(pop) %>% count()

ggplot(pca_cocl, aes(x = PC1,y = PC2,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC1 (3.81%)")+
  ylab("PC2 (1.19%)")+
  scale_color_manual(values = c("#C62828FF",
                                "#F44336FF",
                                "#673AB7FF",
                                "#3F51B5FF",
                                "#2196F3FF",
                                "#4CAF50FF",
                                "#FFEB3BFF",
                                "#FF9800FF",
                                "#9E9E9EFF"), 
                     name = "Population (n)",
                     label = c("KAA",
                               "ROG",
                               "LAG",
                               "BEA",
                               "SAB",
                               "EAS",
                               "RUP",
                               "NOT",
                               "Unassigned"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("figure/reassigned_pc1_pc2_494k_snp_ld_pruned.jpeg",height=8,width = 10)

ggplot(pca_cocl, aes(x = PC3,y = PC4,color = fct_reorder(pop,desc(lat))))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC3 (0.96%)")+
  ylab("PC4 (0.82%)")+
  scale_color_manual(values = c("#C62828FF",
                                "#F44336FF",
                                "#673AB7FF",
                                "#3F51B5FF",
                                "#2196F3FF",
                                "#4CAF50FF",
                                "#FFEB3BFF",
                                "#FF9800FF",
                                "#9E9E9EFF"), 
                     name = "Population (n)",
                     label = c("KAA",
                               "ROG",
                               "LAG",
                               "BEA",
                               "SAB",
                               "EAS",
                               "RUP",
                               "NOT",
                               "Unassigned"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("figure/reassigned_pc3_pc4_494k_snp_ld_pruned.jpeg",height=8,width = 10)

#####################################
#LDA analysis


deme_df<-pca_cocl %>% filter(!is.na(pop))

lda1 <- lda(pop ~ PC1+PC2+PC3+PC4,data=deme_df)

deme_values <- predict(lda1)

lda_values <- as.data.frame(deme_values$x)

lda_values$pop <- pca_cocl %>% filter(!is.na(pop)) %>% pull(pop)

lda_values$lat <- pca_cocl %>% filter(!is.na(pop)) %>% pull(lat)

lda_values %<>% pivot_longer(LD1:LD4,names_to = "LD",values_to = "disc_value")

lda_values %>% 
  ggplot(aes(x=disc_value,y=fct_reorder(pop,lat),fill= fct_reorder(pop,desc(lat))))+
  geom_density_ridges(alpha = 0.6)+
  xlab("Discriminant function")+
  ylab("Deme")+
  scale_fill_manual(values = c("#C62828FF",
                                "#F44336FF",
                                "#673AB7FF",
                                "#3F51B5FF",
                                "#2196F3FF",
                                "#4CAF50FF",
                                "#FFEB3BFF",
                                "#FF9800FF",
                                "#9E9E9EFF"))+
  facet_grid(~LD,scales = "free_x")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


###########################
#####LDA on all dataset####
###########################

deme_df<-pca_cocl %>% filter(!is.na(pop))

deme_df$pop <- as.factor(deme_df$pop)

############################
#####Prediction accuracy####
############################

results <-  matrix(ncol=8, nrow = 100)
reslist <- vector("list", 3)

prop_sample<-c(0.5,0.75,0.95)

for(i in 1:length(prop_sample)){
  for(j in 1:100){
    #Split the data into training (80%) and test set (20%)
    training.individuals <- deme_df$pop %>% 
      createDataPartition(p = prop_sample[i], list = FALSE)
    train.data <- deme_df[training.individuals, ]
    test.data <- deme_df[-training.individuals, ]
    
    # Estimate preprocessing parameters
    preproc.parameter <- train.data %>% 
      preProcess(method = c("center", "scale"))
    
    # Transform the data using the estimated parameters
    train.transform <- preproc.parameter %>% predict(train.data)
    test.transform <- preproc.parameter %>% predict(test.data)
    
    # Fit the model
    model <- lda(pop ~ PC1+PC2+PC3+PC4, data = train.transform)
    
    # Make predictions
    predictions <- model %>% predict(test.transform)
    
    conf <- table(list(predicted=predictions$class, observed=test.data$pop))
    
    results[j,] = diag(conf) / rowSums(conf)
  }
  reslist[[i]] <- results
}


df_res<-as.data.frame(do.call(rbind.data.frame,  reslist))

df_res$prop <- c(rep("0.5",100),rep("0.75",100),rep("0.95",100))

colnames(df_res) <- c("BEA","EAS","KAA","LAG","NOT","ROG","RUP","SAB","prop")

df_res %<>% pivot_longer(BEA:SAB,names_to = "pop",values_to="assignment_acc")

df_res %<>% left_join(pop_file,by=c("pop"="pop_code"))

df_res %>% 
  ggplot(aes(x=fct_reorder(pop,desc(lat)),y=assignment_acc,fill=prop))+
  geom_boxplot()+
  scale_fill_manual(values = c("white","lightgrey","darkgrey"),
                    name="Proportion of\nsamples in training set")+
  xlab("")+
  ylab("Assignment accurary")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18),
        legend.title=element_text(size = 18),
        plot.title= element_text(size = 18))

ggsave("figure/prop_assigned_reassigned.jpeg",width=15,height=10)








