rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(ggh4x)
library(gghalves)
library(patchwork)
library(multcompView)
library(tukeydar)

#---> function for wilcoxon-test
tri.to.squ<-function(x)
{
  rn<-row.names(x)
  cn<-colnames(x)
  an<-unique(c(cn,rn))
  myval<-x[!is.na(x)]
  mymat<-matrix(1,nrow=length(an),ncol=length(an),dimnames=list(an,an))
  for(ext in 1:length(cn))
  {
    for(int in 1:length(rn))
    {
      if(is.na(x[row.names(x)==rn[int],colnames(x)==cn[ext]])) next
      mymat[row.names(mymat)==rn[int],colnames(mymat)==cn[ext]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
      mymat[row.names(mymat)==cn[ext],colnames(mymat)==rn[int]]<-x[row.names(x)==rn[int],colnames(x)==cn[ext]]
    }
    
  }
  return(mymat)
}


#---------------#
#----- Data ----#
#---------------#

het_cocl <- read.table("data/heterozygosity/het_cocl.tsv",fill = T)

whitefish_assignment <- read.table("data/dispersal_model/whitefish_assignment.txt",header=T)

whitefish_assignment %<>% mutate(site = case_when(str_detect(id,"ROG") ~ "ROG",
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

whitefish_assignment %<>% mutate(best_lik2 = case_when(str_detect(best_lik,"Roggan") ~ "ROG",
                                                       str_detect(best_lik,"Kaapsaoui") ~ "KAA",
                                                       str_detect(best_lik,"Biagadwshi") ~ "BIA-L",
                                                       str_detect(best_lik,"La_Grande") ~ "BIA-L",
                                                       str_detect(best_lik,"Beaver") ~ "BEA",
                                                       str_detect(best_lik,"Sabacunica") ~ "SAB",
                                                       str_detect(best_lik,"Eastmain") ~ "EAS-C",
                                                       str_detect(best_lik,"Rupert") ~ "RUP",
                                                       str_detect(best_lik,"Nottaway") ~ "NOT"))

whitefish_assignment %<>% mutate(id2 = gsub("cocl","",whitefish_assignment$id))

het_cocl %<>% left_join(whitefish_assignment,by = c("V1"="id2"))

cocl_f <- read.table("data/heterozygosity/ngsf_cocl")

het_cocl$f <- cocl_f$V1


het_safo <- read.table("data/heterozygosity/het_safo.tsv")


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

trout_assignment %<>% mutate(best_lik2 = case_when(str_detect(best_lik,"Roggan") ~ "ROG-K",
                                                   str_detect(best_lik,"Kaapsaoui") ~ "ROG-K",
                                                   str_detect(best_lik,"Biagadwshi") ~ "BIA",
                                                   str_detect(best_lik,"La_Grande") ~ "LAG",
                                                   str_detect(best_lik,"Aquatuc") ~ "AQU",
                                                   str_detect(best_lik,"Beaver") ~ "BEA",
                                                   str_detect(best_lik,"Sabacunica") ~ "SAB",
                                                   str_detect(best_lik,"Renoyer") ~ "REN",
                                                   str_detect(best_lik,"Conn") ~ "CON",
                                                   str_detect(best_lik,"Eastmain") ~ "EAS",
                                                   str_detect(best_lik,"Broadback") ~ "BRO"))


trout_assignment %<>% mutate(id2 = gsub("safo","",trout_assignment$id))

het_safo %<>% left_join(trout_assignment,by = c("V1"="id2"))

safo_f <- read.table("data/heterozygosity/ngsf_safo")

het_safo$f <- safo_f$V1

het <- rbind(het_cocl,het_safo)

het$sp <- c(rep("cocl",456),rep("safo",362))

het$site <- factor(het$site,levels = c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","JAC","RUP","BRO","NOT"))



#---> Heterozygosity


#Whitefish

whitefish <- ggplot(het %>% filter(!is.na(site),sp=="cocl"),aes(x = site, y = V2,fill=site))+
  geom_boxplot(
    width = .4, position= position_nudge(x=.2),
    outlier.shape = NA)+
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5,
    aes(color = site),
    size = 2
  ) +
  ylim(c(0.0015,0.0045))+
  xlab("Location")+
  ylab("")+
  scale_fill_manual(values = rev(c("#A50026FF",
                                   "#DD3D2DFF",
                                   "#fddbc7",
                                   "#F67E4BFF",
                                   "#FDB366FF",
                                   "#FEDA8BFF",
                                   "#EAECCCFF",
                                   "#abdda4",
                                   "#74add1",
                                   "#4575b4",
                                   "#313695",
                                   "#5e4fa2")
  ))+
  scale_color_manual(values = rev(c("#A50026FF",
                                    "#DD3D2DFF",
                                    "#fddbc7",
                                    "#F67E4BFF",
                                    "#FDB366FF",
                                    "#FEDA8BFF",
                                    "#EAECCCFF",
                                    "#abdda4",
                                    "#74add1",
                                    "#4575b4",
                                    "#313695",
                                    "#5e4fa2")
  ))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle=45,vjust = .9,hjust=1),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))

whitefish

ggsave("figures/het_whitefish.png")

het$site2 <- factor(het$site,levels = rev(c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","JAC","RUP","BRO","NOT")))

het$best_lik2 <- factor(het$best_lik2,levels = rev(c("ROG","KAA", "ROG-K","BIA-L","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","EAS-C","RUP","BRO","NOT")))

whitefish2 <- ggplot(het %>% filter(!is.na(best_lik2),sp=="cocl"),aes(x = best_lik2, y = V2))+
  geom_boxplot(
    width = .4, position= position_nudge(x=.2),
    outlier.shape = NA)+
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5,
    size = 2
  ) +
  ylim(c(0.0015,0.0045))+
  ylab("")+
  xlab("")+
  coord_flip()+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))

whitefish2

#---> Trout

trout <- ggplot(het %>% filter(!is.na(site),sp=="safo"),aes(x = site, y = V2,fill=site))+
  geom_boxplot(
    width = .4, position= position_nudge(x=.2),
    outlier.shape = NA)+
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5,
    aes(color = site),
    size = 2
  ) +
  ylim(c(0.0015,0.0045))+
  xlab("Location")+
  ylab("")+
  scale_fill_manual(values = rev(c(
    "#d53e4f",
    "#F67E4BFF",
    "#FDB366FF",
    "#FEDA8BFF",
    "#EAECCCFF",
    "#abdda4",
    "#66c2a5",
    "#74add1",
    "#4575b4",
    "#313695",
    "#5e4fa2")
  ))+
  scale_color_manual(values = rev(c(
    "#d53e4f",
    "#F67E4BFF",
    "#FDB366FF",
    "#FEDA8BFF",
    "#EAECCCFF",
    "#abdda4",
    "#66c2a5",
    "#74add1",
    "#4575b4",
    "#313695",
    "#5e4fa2")
  ))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle=45,vjust = .9,hjust=1),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))

trout

ggsave("figures/het_trout.png")

trout2 <- ggplot(het %>% filter(!is.na(best_lik2),sp=="safo"),aes(x = best_lik2, y = V2))+
  geom_boxplot(
    width = .4, position= position_nudge(x=.2),
    outlier.shape = NA)+
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5,
    size = 2
  ) +
  xlab("")+
  ylab("Multi-locus heterozygosity")+
  ylim(c(0.0015,0.0045))+
  coord_flip()+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))


trout2

#Species comp

sp_comp <- ggplot(het %>% filter(!is.na(site)),aes(x = sp, y = V2))+
  geom_boxplot(width = .2, 
               outlier.shape = NA, notch = T)+
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .2, 
    alpha = .2,
    size = 2
  ) +
  scale_x_discrete(labels=c("cocl" = expression(italic("C. clupeaformis")), 
                            "safo" = expression(italic("S. fontinalis"))))+
  ylim(c(0.0015,0.007))+
  xlab("")+
  ylab("Multi-locus heterozygosity")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))

sp_comp

ggsave("figures/species_het.png")

ggplot(het %>% filter(!is.na(site)),aes(x = sp, y = asinh(f)))+
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .2, 
    alpha = .2,
    size = 2
  ) +
  scale_x_discrete(labels=c("cocl" = expression(italic("C. clupeaformis")), 
                            "safo" = expression(italic("S. fontinalis"))))+
  xlab("")+
  ylab("Inbreeding coefficient (F)")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))


qqplot(safo_f$V1,cocl_f$V1)
abline(c(0,1))

eda_qq(safo_f$V1,cocl_f$V1)


wilcox.test(het$f ~ het$sp)


############################
#---> Kruskal-Wallis test 

# Whitefish 

het_cocl_aov <- 
  het %<>% 
  filter(!is.na(site),sp=="cocl") %>%
  dplyr::select(V1,V2,site) %>% 
  dplyr::rename(id = V1,
         het = V2)


het_cocl_aov %>% group_by(site) %>% dplyr::summarise(mean_het = mean(het, na.rm =T))


kw_cocl <- kruskal.test(het ~ site, data = het_cocl_aov)

#Kruskal-Wallis chi-squared = 136.62, df = 11, p-value < 2.2e-16

wt_cocl <- pairwise.wilcox.test(het_cocl_aov$het, het_cocl_aov$site, paired = F,
                     p.adjust.method = "BH")


mymat<-tri.to.squ(wt_cocl$p.value)

letters_cocl<-multcompLetters(mymat)
 
group_cocl <- as_tibble(as.list(letters_cocl$Letters)) %>% pivot_longer(cols = ROG:NOT,names_to = "site",values_to = "group_wt")

#Figure with groups

cocl <- het %>% filter(!is.na(site),sp=="cocl")

cocl %<>% left_join(group_cocl,by = "site")

cocl %<>% mutate(site2 = case_when(is.na(site) & str_detect(V1,"BEA") ~ "BEA",
                                  is.na(site) & str_detect(V1,"SAB") ~ "SAB",
                                  is.na(site) & str_detect(V1, "LAG") ~ "LAG",
                                  T ~ site))

cocl$site2 <- factor(cocl$site2,levels = c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","JAC","RUP","BRO","NOT"))


whitefish <- ggplot(cocl,aes(x = site2, y = V2,fill=site2))+
  geom_boxplot(
    width = .4, position= position_nudge(x=.2),
    outlier.shape = NA)+
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5,
    aes(color = site2),
    size = 2
  ) +
  geom_text(data=cocl %>% group_by(site2) %>% summarise(y = min(max(V2), diff(quantile(V2, probs = c(0.25, 0.75),na.rm=T) * c(3/2, 5/2)),na.rm=T),group_wts = first(group_wt)),
            aes(y=y,label = group_wts),
            vjust = -0.5,
            size = 5)+
  ylim(c(0.0015,0.0045))+
  xlab("Location")+
  ylab("")+
  scale_fill_manual(values = rev(c("#A50026FF",
                                   "#DD3D2DFF",
                                   "#fddbc7",
                                   "#F67E4BFF",
                                   "#FDB366FF",
                                   "#FEDA8BFF",
                                   "#EAECCCFF",
                                   "#abdda4",
                                   "#74add1",
                                   "#4575b4",
                                   "#313695",
                                   "#5e4fa2")
  ))+
  scale_color_manual(values = rev(c("#A50026FF",
                                    "#DD3D2DFF",
                                    "#fddbc7",
                                    "#F67E4BFF",
                                    "#FDB366FF",
                                    "#FEDA8BFF",
                                    "#EAECCCFF",
                                    "#abdda4",
                                    "#74add1",
                                    "#4575b4",
                                    "#313695",
                                    "#5e4fa2")
  ))+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle=45,vjust = .9,hjust=1),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))


ggsave("figures/het_whitefish.png")


# Trout

het_safo_aov <- 
  het %>% 
  filter(!is.na(site),sp=="safo") %>%
  dplyr::select(V1,V2,site) %>% 
  dplyr::rename(id = V1,
                het = V2)

het_safo_aov %>% group_by(site) %>% dplyr::summarise(mean_het = mean(het, na.rm =T))


kw_safo <- kruskal.test(het ~ site, data = het_safo_aov)

#Kruskal-Wallis chi-squared = 119.81, df = 10, p-value < 2.2e-16

wt_safo <- pairwise.wilcox.test(het_safo_aov$het, het_safo_aov$site, paired = F,
                                p.adjust.method = "BH")


mymat<-tri.to.squ(wt_safo$p.value)

letters_safo<-multcompLetters(mymat)

group_wt_safo <- as_tibble(as.list(letters_safo$Letters)) %>% pivot_longer(cols = ROG:BRO,names_to = "site",values_to = "group_wt")

#Figure with groups

safo <- het %>% filter(!is.na(site),sp=="safo")

safo %<>% left_join(group_wt_safo,by = "site")

safo$site <- factor(safo$site,levels = c("ROG","KAA","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","BRO"))

trout <- ggplot(safo,aes(x = site, y = V2,fill=site))+
  geom_boxplot(
    width = .4, position= position_nudge(x=.2),
    outlier.shape = NA)+
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5,
    aes(color = site),
    size = 2
  ) +
  geom_text(data=safo %>% group_by(site) %>% summarise(y = min(max(V2), diff(quantile(V2, probs = c(0.25, 0.75),na.rm=T) * c(3/2, 5/2)),na.rm=T),group_wts = first(group_wt)),
            aes(y=y,label = group_wts),
            vjust = -0.5,
            size = 5)+
  scale_fill_manual(values = rev(c(
    "#d53e4f",
    "#F67E4BFF",
    "#FDB366FF",
    "#FEDA8BFF",
    "#EAECCCFF",
    "#abdda4",
    "#66c2a5",
    "#74add1",
    "#4575b4",
    "#313695",
    "#5e4fa2")
  ))+
  scale_color_manual(values = rev(c(
    "#d53e4f",
    "#F67E4BFF",
    "#FDB366FF",
    "#FEDA8BFF",
    "#EAECCCFF",
    "#abdda4",
    "#66c2a5",
    "#74add1",
    "#4575b4",
    "#313695",
    "#5e4fa2")
  ))+
  ylim(c(0.0015,0.0045))+
  xlab("Location")+
  ylab("")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18,angle=45,vjust = .9,hjust=1),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "none",
        strip.text = element_text(size = 18),
        strip.background = element_rect(fill = "white"),panel.spacing = unit(2, "lines"))



#--->

sp_comp + whitefish + trout + 
  plot_layout(widths = c(0.4,0.6,0.6))

ggsave("figures/het_spcomp_cocl_safo.png",width = 20,height = 7)
ggsave("figures/het_spcomp_cocl_safo.pdf",width = 20,height = 7)


t.test(het$V2~het$sp)




