rm(list = ls())

#####################
#######Library#######
#####################

library(tidyverse)
library(magrittr)
library(paletteer)

#####################
#########DATA########
#####################

samples <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/ind_list_analysis/ind.list")

colnames(samples) <- "id"

samples %<>% mutate(pop = case_when(str_detect(string = id, "RUP") ~ "RUP",
                                     str_detect(string = id, "NOT") ~ "NOT",
                                     str_detect(string = id, "BEA") ~ "BEA",
                                     str_detect(string = id, "FIS") ~ "EAS",
                                     str_detect(string =id, "JAC") ~ "JAC",
                                     str_detect(string = id, "SAB") ~ "SAB",
                                     str_detect(string = id, "BIA") ~ "BIA",
                                     str_detect(string = id, "LAG") ~ "LAG",
                                     str_detect(string = id, "TIL") ~ "REN",
                                     str_detect(string = id, "MOA") ~ "MOA",
                                     str_detect(string = id, "ROG") ~ "ROG",
                                     str_detect(string = id, "EAS") ~ "EAS",
                                     str_detect(string = id, "CHIm") ~ "KAA",
                                     str_detect(string = id, "CON") ~ "CON",
                                     str_detect(string = id, "BRO") ~ "BRO",
                                     str_detect(string = id, "WAS") ~ "WAS",
                                     str_detect(string = id, "REN") ~ "REN",
                                     str_detect(string = id, "AQU") ~ "AQU"))

pop_file <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/pop_sampled.txt",header=T)

samples %<>% left_join(pop_file,by = c("pop" = "pop_code"))

missing <- read.table("data/missing_data/missing_per_ind.txt")

samples$missing <- missing$V2


ggplot(samples, aes(x = fct_reorder(pop,desc(lat)), y = missing))+
  geom_boxplot()+
  xlab("Population (ordered by latitude)")+
  ylab("Proportion of missing data")+
  theme_bw()

ggsave("figure/missing_data_per_pop.jpeg")

mean(samples$missing)
mean(samples$missing) + 1.96*sd(samples$missing)


ggplot(samples,aes(x=missing))+
  geom_histogram(color = "black",bins = 60)+
  geom_vline(xintercept = mean(samples$missing) + 2*sd(samples$missing))+
  theme_bw()

samples[which(samples$missing > mean(samples$missing) + 2*sd(samples$missing)),]

samples_to_keep <- samples %>% filter(missing < mean(samples$missing) + 2*sd(samples$missing)) %>% pull(id)

samples_to_keep_final <- unique(samples_to_keep)

write.table(samples_to_keep_final,"data/ind_list_without_missing.txt",quote = F,row.names = F,col.names = F)






