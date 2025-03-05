rm(list = ls())

#####################
#######Library#######
#####################

library(tidyverse)
library(magrittr)
library(paletteer)
library(plotly)
library(qvalue)
library(qqman)
library(bigutilsr)
library(RcppCNPy)
library(data.table)

#####################
#########DATA########
#####################

setwd("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/pcadapt/")

z <- npyLoad("all_maf0.05pctind0.85_maxdepth10_ALL_CHR_singletons.pcadapt.pcadapt.zscores.npy")
sites <- fread("sites_singletons_all_chr")

# For one component only
d <- dist_ogk(z[,1:4])

p <- pchisq(d,3,lower.tail = F)

pcadapt <- bind_cols(sites,p)

#Change col names
colnames(pcadapt) <- c("chr","pos","maj","min","pvalue1")

outliers <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/data/gea/gea_canditates_unconditioned_rda_all_snp.tsv",header = T)

#########################
#Format for manhattan plot

pcadapt %<>% mutate(chr_num = as.numeric(str_replace(chr,"Chr","")))

pcadapt %<>% mutate(ID = row_number(pcadapt))

outliers %<>% separate(snp,into = c("chr","pos"))

outliers$pos <- as.numeric(outliers$pos)

outliers_id <- outliers %>% left_join(pcadapt, by = c("chr"="chr","pos"="pos")) %>% pull(ID)

png("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/figure/pcadapt_manhattant_pc2_pc4.jpeg",height=30,width = 50,units = "cm",res=300)

manhattan(pcadapt, chr="chr_num", bp="pos", snp="ID", p="pvalue1")

dev.off()

png("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/figure/pcadapt_manhattan_GEA_outliers.jpeg",height=30,width = 50,units = "cm",res=300)

manhattan(pcadapt, chr="chr_num", bp="pos", snp="ID", p="pvalue1", highlight = outliers_id)

dev.off()

########################

outliers_pcadapt <- pcadapt %>% filter(-log10(pvalue1) > 7) %>% dplyr::select(chr,pos,pvalue1)

write.table(outliers_pcadapt, "pcadapt_outliers.txt",row.names = F,col.names = F,quote = F)

woutliers_pcadapt_chr35 <- pcadapt %>% filter(-log10(pvalue1) > 20, chr == "Chr35") %>% select(chr,pos,pvalue1)

write.table(outliers_pcadapt_chr35, "data/pcadapt/pcadapt_outliers_chr35.txt",row.names = F,col.names = F,quote = F)

outliers_pcadapt_chr7 <- pcadapt %>% filter(-log10(pvalue1) > 7, pos < 122600000, chr == "Chr07") %>% select(chr,pos,pvalue1)

write.table(outliers_pcadapt_chr7, "data/pcadapt/pcadapt_outliers_chr7.txt",row.names = F,col.names = F,quote = F)


overlap <- outliers_pcadapt %>% left_join(outliers,by = c("chr","pos"))

overlap %>% filter(!is.na(loading))

#

outliers <- read.table("data/pcadapt/pcadapt_outliers.txt")

out_chr <- outliers %>% group_by(V1) %>% count() 

out_chr %<>% mutate(perc = n/sum(out_chr$n)*100) %>% arrange(desc(perc))

out_chr %>% ggplot(aes(x = reorder(V1,desc(perc)), y = perc))+
  geom_bar(stat = "identity")+
  ylab("Proportion of outliers in chromosome (%)")+
  xlab("")+
  theme_bw()+
  theme(
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18, angle = 90),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none")

ggsave("figure/prop_outlier_per_chr.jpeg",width = 12, height = 8)















