rm(list = ls())

library(tidyverse)
library(magrittr)
library(paletteer)

pca <- read.table("~/Desktop/Baie_James_Paper/whitefish_pop_gen/04_pca_cerri/all_maf0.05_pctind0.80_maxdepth10.cov.pca")

sample <- read.table("04_pca_cerri/bam.filelist")

pca$sample_id <- sample$V1

pca %<>% mutate(pop = case_when(str_detect(string = sample_id, "RUP") ~ "RUP",
                                 str_detect(string = sample_id, "NOT") ~ "NOT",
                                 str_detect(string = sample_id, "BEA") ~ "BEA",
                                 str_detect(string = sample_id, "FIS") ~ "FIS",
                                 str_detect(string = sample_id, "JAC") ~ "JAC",
                                 str_detect(string = sample_id, "SAB") ~ "SAB",
                                 str_detect(string = sample_id, "BIA") ~ "BIA",
                                 str_detect(string = sample_id, "LAG") ~ "LAG",
                                 str_detect(string = sample_id, "TIL") ~ "TIL",
                                 str_detect(string = sample_id, "MOA") ~ "SAB",
                                str_detect(string = sample_id, "ROG") ~ "ROG"))

eig <- read.table("~/Desktop/Baie_James_Paper/whitefish_pop_gen/04_pca_cerri/all_maf0.05_pctind0.80_maxdepth10.cov.eig")

head(eig)

ggplot(pca, aes(x = PC1,y = PC2,color = pop))+
  geom_point(alpha=0.8,size = 3)+
  xlab("PC1 (13.94%)")+
  ylab("PC2 (1.75%)")+
  scale_color_manual(values = paletteer_d(`"awtools::bpalette"`))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

cisco <- pca %>% filter(PC1 < -0.5) %>% pull(sample_id)

write.table(cisco, "cisco_samples.txt",quote=F,row.names = F,col.names = F)


