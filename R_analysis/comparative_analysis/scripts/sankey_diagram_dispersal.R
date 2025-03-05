rm(list=ls())

setwd(dir = "~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/")

#---------------#
#--- Library ---#
#---------------#

library(tidyverse)
library(magrittr)
library(patchwork)
library(plyr)
library(ggsankey)
library(networkD3)
library(betareg)
library(lmtest)
library(car)
library(effects)
library(jtools)
library(performance)


#---------------#
#----- Data ----#
#---------------#

#---> TROUT

trout_assignment <- read.table("data/dispersal_model/trout_assignment.txt",header=T)

trout_assignment %<>% mutate(site = case_when(str_detect(id,"ROG") ~ "ROG-K",
                                             str_detect(id,"CHIm") ~ "ROG-K",
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

trout_assignment %<>% mutate(site = if_else(str_detect(id,"WAS"),best_lik2,site))


#Create a column which identifies individual as migrant

trout_assignment %<>% mutate(migrant = if_else(best_lik2 == site,"no","yes"))

trout_assignment %>% group_by(migrant) %>% dplyr::count()

trout_assignment %>% filter(post_prob >= 0.8) %>% group_by(migrant) %>% dplyr::count()

trout_assignment %<>% filter(post_prob >= 0.8)

#Sankey diagram

trout_sankey <- as.data.frame.matrix(table(trout_assignment$best_lik2,trout_assignment$site))

# I need a long format
data_long <- trout_sankey %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)

colnames(data_long) <- c("source", "target", "value")

data_long$target <- paste(data_long$target, " ", sep="")

data_long$source <- factor(data_long$source,levels = c("ROG-K","BIA","LAG","AQU","BEA","SAB","REN","CON","EAS","BRO"))
data_long$target <- factor(data_long$target,levels = c("ROG-K ","BIA ","LAG ","AQU ","BEA ","SAB ","REN ","CON ","EAS ","BRO "))


# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name=c(levels(data_long$source), levels(data_long$target)) %>% unique())

nodes$group <- rep(c(levels(data_long$source)),2)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long$IDsource=match(data_long$source, nodes$name)-1 
data_long$IDtarget=match(data_long$target, nodes$name)-1


# prepare colour scale

ColourScal ='d3.scaleOrdinal() .range(["#66c2a5", "#abdda4","#4575b4","#313695", "#d53e4f", "#FDB366FF", "#F67E4BFF", "#FEDA8BFF", "#74add1","#EAECCCFF"])'


# Make the Network

sankeyNetwork(Links = data_long, 
              Nodes = nodes,
              Source = "IDsource", 
              Target = "IDtarget",
              Value = "value", 
              NodeID = "name", 
              sinksRight=FALSE, 
              colourScale=ColourScal,
              LinkGroup = "source",
              NodeGroup="group",
              nodeWidth=10, 
              fontSize=13, 
              nodePadding=20,fontFamily = "arial")


#---> Whitefish

whitefish_assignment <- read.table("data/dispersal_model/whitefish_assignment.txt",header=T)

whitefish_assignment %<>% mutate(site = case_when(str_detect(id,"ROG") ~ "ROG",
                                              str_detect(id,"CHIm") ~ "KAA",
                                              str_detect(id,"BIA") ~ "BIA-L",
                                              str_detect(id,"LAG") & best_lik == "La_Grande" ~ "BIA-L",
                                              str_detect(id,"LAG") & best_lik == "Roggan" ~ "ROG",
                                              str_detect(id,"BEA") ~ "BEA",
                                              str_detect(id,"SAB") ~ "SAB",
                                              str_detect(id,"MOA") ~ "SAB",
                                              str_detect(id,"REN") ~ "EAS-C",
                                              str_detect(id,"TIL") ~ "EAS-C",
                                              str_detect(id,"CON") ~ "EAS-C",
                                              str_detect(id,"EAS") ~ "EAS-C",
                                              str_detect(id,"FIS") ~ "EAS-C",
                                              str_detect(id,"JAC") ~ "EAS-C",
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

#Create a column which identifies individual as migrant

whitefish_assignment %<>% mutate(migrant = if_else(best_lik2 == site,"no","yes")) 

whitefish_assignment %>% group_by(migrant) %>% dplyr::count()

whitefish_assignment %>% filter(post_prob >= 0.8) %>% group_by(migrant) %>% dplyr::count()

whitefish_assignment %<>% filter(post_prob >= 0.8)

#Sankey diagram

whitefish_sankey <- as.data.frame.matrix(table(whitefish_assignment$best_lik2,whitefish_assignment$site))

# I need a long format
data_long2 <- whitefish_sankey %>%
  rownames_to_column %>%
  gather(key = 'key', value = 'value', -rowname) %>%
  filter(value > 0)

colnames(data_long2) <- c("source", "target", "value")

data_long2$target <- paste(data_long2$target, " ", sep="")

data_long2$source <- factor(data_long2$source,levels = c("ROG","KAA","BIA-L","BEA","SAB","EAS-C","RUP","NOT"))
data_long2$target <- factor(data_long2$target,levels = c("ROG ","KAA ","BIA-L ","BEA ","SAB ","EAS-C ","RUP ","NOT "))


# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes2 <- data.frame(name=c(levels(data_long2$source), levels(data_long2$target)) %>% unique())

nodes2$group <- rep(c(levels(data_long2$source)),2)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long2$IDsource=match(data_long2$source, nodes2$name)-1 
data_long2$IDtarget=match(data_long2$target, nodes2$name)-1


# prepare colour scale

ColourScal ='d3.scaleOrdinal() .range(["#abdda4","#EAECCCFF","#74add1","#fdae61","#a50026","#f46d43","#4575b4","#364B9AFF"])'

"#EAECCCFF","#abdda4","#74add1","#fdae61","#a50026","#f46d43","#4575b4","#364B9AFF"

# Make the Network
sankeyNetwork(Links = data_long2, 
              Nodes = nodes2,
              Source = "IDsource", 
              Target = "IDtarget",
              Value = "value", 
              NodeID = "name", 
              sinksRight=FALSE, 
              colourScale=ColourScal,
              LinkGroup = "source",
              NodeGroup="group",
              nodeWidth=10, 
              fontSize=13, 
              nodePadding=20,fontFamily = "arial")


#Proportion of disperser per population

siteswf <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/pop_dispersal.txt",header=T)


whitefish_disp <- whitefish_assignment %>% 
  mutate(site2 = case_when(str_detect(id,"ROG") ~ "ROG",
                          str_detect(id,"CHIm") ~ "KAA",
                          str_detect(id,"BIA") ~ "BIA",
                          str_detect(id,"LAG") & best_lik == "La_Grande" ~ "LAG",
                          str_detect(id,"LAG") & best_lik == "Roggan" ~ "ROG",
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
                          str_detect(id,"NOT") ~ "NOT")) %>% 
  group_by(site2) %>%
  filter(post_prob >= 0.8) %>% 
  dplyr::count(migrant) %>% 
  dplyr::summarize(
    total_yes = sum(n[migrant == "yes"]),
    total = sum(n),
    proportion_yes = total_yes / total
  )

wf_disp_lat <- whitefish_disp %>% left_join(siteswf, by = c("site2" = "pop_code"))


ggplot(wf_disp_lat, aes(x = lat, y = proportion_yes))+
  geom_point()+
  geom_smooth(method = "lm",color = "black")+
  theme_bw()



trout_disp <- trout_assignment %>%  mutate(site2 = case_when(str_detect(id,"ROG") ~ "ROG",
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
                                                str_detect(id,"WAS") ~ "JAC",
                                                str_detect(id,"BRO") ~ "BRO")) %>% 
  group_by(site2) %>%
  filter(post_prob >= 0.8) %>% 
  dplyr::count(migrant) %>% 
  dplyr::summarize(
    total_yes = sum(n[migrant == "yes"]),
    total = sum(n),
    proportion_yes = total_yes / total
  )

bt_disp_lat <- trout_disp %>% left_join(siteswf, by = c("site2" = "pop_code"))


ggplot(bt_disp_lat, aes(x = lat, y = proportion_yes))+
  geom_point()+
  geom_smooth(method = "lm",color = "black")+
  theme_bw()



disp_lat <- bind_rows(bt_disp_lat,wf_disp_lat) %>% filter(!is.na(site2))

disp_lat$sp <- c(rep("safo",12),rep("cocl",12))

ggplot(disp_lat, aes(x = lat, y = proportion_yes,shape=sp,group=1))+
  geom_point(size=3)+
  geom_smooth(method="glm",se = F,color = "black",linetype=2)+
  scale_shape_manual(values = c(1,16))+
  xlab("Latitude")+
  ylab("Proportion of dispersers")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

disp_lat %>% group_by(sp) %>% dplyr::summarize(mean_im = mean(proportion_yes),
                                               sd_im = sd(proportion_yes))

########################
# Environmental data all increase vif

env_wf <- read.table("../environmental_data/results/whitefish_env_data.txt", header = T)
env_bt <- read.table("../environmental_data/results/trout_env_data.txt", header = T)

env_bt %<>% filter(!SITE_CODE %in% c("MOA","WAS"))

env <- bind_rows(env_wf,env_bt) %>% group_by(SITE_CODE) %>% slice(1)

salinity <- read.table("data/salinity_per_site_Raph.txt", header = T)

env_sal <- left_join(env, salinity, by = c("SITE_CODE" = "pop_code"))

########################
# Test if proportion of dispersers varies with latitude

disp_lat %<>% mutate(local = total-total_yes)

trials = cbind(disp_lat$total_yes, disp_lat$local)


model.log = glm(trials ~ lat + sp ,
                data = disp_lat,
                family = binomial(link="logit"),
                weights = n)


summary(model.log)

exp(coef(model.log))

summ(model.log)

#Check model assumption using boxtidwell
boxTidwell(formula = trials ~ lat,
           data = disp_lat)

#All is good, the term is non significant

# Perform Hosmer-Lemeshow test

performance_hosmer(model.log, n_bins = 3)

model.lognull = glm(trials ~ 1,
                data = disp_lat,
                family = binomial(link="logit"))

summary(model.lognull)

lrtest(model.lognull,model.log,model.log2)

with(summary(model.log), 1 - deviance/null.deviance)

#McFadden’s R-squared = 39%

plot(residuals(model.log , type = "deviance"))

model.log2 = glm(trials ~ lat,
                data = disp_lat,
                family = binomial(link="logit"))

# Get logit predictions of whitefish
logit.predictions <- predict(object = model.log2, se=T)

# Apply inverse logit to transform to probabilities
prob.predictions <- plogis(logit.predictions$fit)

LL <- plogis(logit.predictions$fit - 1.96*logit.predictions$se.fit)
UL <- plogis(logit.predictions$fit + 1.96*logit.predictions$se.fit)

confint()

disp_lat$pred <- prob.predictions
disp_lat$LL <- LL
disp_lat$UL <- UL

ggplot(disp_lat)+
  geom_point(aes(x = lat, y = proportion_yes,shape=sp),size=3)+
  geom_line(aes(x = lat, y = pred), color = "black")+
  geom_ribbon(aes(x = lat, ymin = LL, ymax = UL), alpha = 0.2)+ 
  scale_shape_manual(values = c(16,1),labels = c("",""))+
  xlab("Latitude")+
  ylab("Proportion of dispersers")+
  theme_bw()+
  theme(axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = c(0.7,0.9),
        legend.background = element_blank())


ggsave("figures/disp_lat.jpeg",width = 8, height = 8)


