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

pcadapt <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/pcadapt/pcadapt.pval.txt")
pos <-read.table("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/pcadapt/singletons_all_chr.txt")

pcadapt$chr <- pos$V1
pcadapt$pos <- pos$V2

#Change col names
colnames(pcadapt) <- c("pvalue","chr","pos")

#Add col with id so we can plot position as a continuous variable
pcadapt %<>% rowid_to_column("ID")

#Create col with adjusted pvalue
pcadapt$padj <- p.adjust(pcadapt$pvalue,method="bonferroni")
alpha <- 0.01 #FDR lower than 1%
outliers <- which(pcadapt$padj < alpha)
length(outliers)

#1125 SNPs identified as outlier with FDR < 1%

#############
#Plot results


pcadapt %<>% mutate(chr_num = as.numeric(str_replace(chr,"LG","")))

pcadapt$chr_col <- pcadapt$chr_num %% 2 == 0 

#########################
#graph avec qqman package

png("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/figures/manhattan_pcadapt.png",height=30,width = 50,units = "cm",res=300)

manhattan(pcadapt, chr="chr_num", bp="pos", snp="ID", p="pvalue")

dev.off()

########################

outliers_pcadapt <- pcadapt %>% filter(-log10(pvalue) > 7) %>% dplyr::select(chr,pos,pvalue)

write.table(outliers_pcadapt, "data/pcadapt/pcadapt_outliers_trout.txt",row.names = F,col.names = F,quote = F)


#

outliers <- read.table("data/pcadapt/pcadapt_outliers_trout.txt")

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

ggsave("figures/prop_outlier_per_chr.jpeg",width = 11, height = 8)



