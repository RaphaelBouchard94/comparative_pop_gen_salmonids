rm(list = ls())

########################
# 1 - Library  #
########################

library(sp)
library(sf)
library(httr)
library(tidyverse)
library(magrittr)
library(qrmtools)
library(ggspatial)
library("rnaturalearth")
library("rnaturalearthdata")
library("maps")
library(wesanderson)
library(patchwork)
library(marmap)
library(jpeg)
library(paletteer)
library(readxl)


########################
# 2 - Load Data  #
########################

# 1. GET RIVERS DATA
#---------

#World map background
world <- map("world", fill=TRUE, plot=FALSE)

world <- ne_countries(
  country = "Canada",
  scale = "large", returnclass = "sf"
)

world_geo <- st_set_crs(world, "EPSG:4326")

#River shapefile data
nam_riv <- st_read("HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na_shp")

nam_riv %<>% 
  mutate(width = as.numeric(ORD_FLOW),
         width = case_when(width == 3 ~ 1,
                           width == 4 ~ 0.8,
                           width == 5 ~ 0.6,
                           width == 6 ~ 0.4,
                           width == 7 ~ 0.2,
                           width == 8 ~ 0.2,
                           width == 9 ~ 0.1,
                           width == 10 ~ 0.1,
                           TRUE ~ 0)) %>%
  st_as_sf()

crop_bbox <- c(xmin = -82.084, ymin = 51.021, xmax = -74.778, ymax = 55)

river_cropped <- st_crop(nam_riv, crop_bbox)

#Bathymetry

bathy <- getNOAA.bathy(-82.084, -79, 51.80, 55, res = 1, keep = TRUE)
ggbathy <- fortify(bathy)


bathy <- getNOAA.bathy(-82.084, -74.778, 51.021, 55, res =1, keep = TRUE)
autoplot(bathy, geom=c("r", "c"), colour="white", size=0.1) + scale_fill_etopo()


#Sampling location data

##Reformatage de long/lat à partir du document original

pop <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/pop_sampled.txt",header = T)

pop

pop_sf <- st_as_sf(pop, coords = c("long", "lat"), crs = 4326)

pop_sf$lat <- pop$lat

#############################

ggplot() +
  geom_sf(data = world_geo, fill = "grey90")+
  geom_sf(data = pop_sf %>% filter(pop != "MoarBay"), 
          aes(fill = fct_reorder(pop,desc(lat))), color = "black",size = 5,pch = 21)+
  coord_sf(c(-81.75, -76), ylim = c(51.25, 54.75))+
  scale_fill_manual(values = paletteer_d(`"awtools::bpalette"`), name = "Population (n)")+
  scale_size(range=c(0, .5)) +
  annotate(geom = "text", 
           x = -80.7, y = 54, 
           label = "James\nBay", 
           color = "lightblue", 
           size = 8,
           fontface = "italic")+
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book"
  ) +
  theme(panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=18),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.text.x = element_text(size = 26),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))

ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/whitefish_pop_gen/figure/map_whitefish.jpeg",width =10, height =10)

##########################
#cisco


pop_cisco <- read.table("~/Desktop/Baie_James_Paper/pop_sampled.txt",header = T)

pop_cisco %<>% filter(pop_code %in% c("LAG","SAB","FIS","RUP","NOT"))

pop_cisco %<>% mutate(pop_code_cisco = if_else(pop_code == "FIS","EAS",as.character(pop_code)))

pop_sf_cisco <- st_as_sf(pop_cisco, coords = c("long", "lat"), crs = 4326)

pop_sf_cisco$lat <- pop_cisco$lat


cols<-c("LAG" = "#9C27B0FF",
        "SAB" = "#2196F3FF",
        "EAS" = "#29B6F6FF",
        "RUP" = "#4CAF50FF",
        "NOT" = "#8BC34AFF")


ggplot() +
  geom_sf(data = world_geo, fill = "grey90")+
  geom_sf(data=river_cropped, aes(size=width),show.legend = F) +
  geom_sf(data = pop_sf_cisco, 
          aes(fill = fct_reorder(pop_code_cisco,desc(lat))), color = "black",size = 5,pch = 21)+
  coord_sf(c(-81.75, -76), ylim = c(51.25, 54.75))+
  scale_fill_manual(labels = c("La Grande",
                               "Sabacunica",
                               "Eastmain",
                               "Rupert",
                               "Nottaway"),
                               values = cols, name = "Population (n)")+
  scale_size(range=c(0, .5)) +
  annotate(geom = "text", 
           x = -80.7, y = 54, 
           label = "James\nBay", 
           color = "blue", 
           size = 8,
           fontface = "italic")+
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book"
  ) +
  theme(panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=18),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.text.x = element_text(size = 26),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("~/Desktop/Baie_James_Paper/map_cisco.jpeg",width =10, height =10)


##########################
#Brook Trout


pop_trout<- read_xlsx("~/Desktop/Baie_James_Paper/01_population_genomics/trout_pop_gen/data/source_pop_bt.xlsx")

pop_sf_trout <- st_as_sf(pop_trout, coords = c("long", "lat"), crs = 4326)

pop_sf_trout$lat <- pop_trout$lat


cols<-c("LAG" = "#9C27B0FF",
        "SAB" = "#2196F3FF",
        "EAS" = "#29B6F6FF",
        "RUP" = "#4CAF50FF",
        "NOT" = "#8BC34AFF")


ggplot() +
  geom_sf(data = world_geo, fill = "grey90")+
  geom_sf(data=river_cropped, aes(size=width),show.legend = F) +
  geom_sf(data = pop_sf_trout, 
          aes(fill = fct_reorder(pop_code,desc(lat))), color = "black",size = 5,pch = 21)+
  coord_sf(c(-81.75, -76), ylim = c(51.25, 54.75))+
  scale_fill_manual(labels = c("Roggan R.(n:16)",
                               "Kaapsaoui R. (n:40) ",
                               "Biagadwshi R. (n:36) ",
                               "La Grande R. (n:39)",
                               "Aquatuc R. (n:21)",
                               "Beaver R. (n:32)",
                               "Sabacunica R. (n:32)",
                               "Renoyer R. (n:32)",
                               "Fishing R. (n:43)",
                               "Broadback R. (n:15)"),
                    values = paletteer_d(`"awtools::bpalette"`), name = "Population (n)")+
  scale_size(range=c(0, .5)) +
  annotate(geom = "text", 
           x = -80.7, y = 54, 
           label = "James\nBay", 
           color = "blue", 
           size = 8,
           fontface = "italic")+
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book"
  ) +
  theme(panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=18),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.text.x = element_text(size = 26),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggplot() +
  geom_sf(data = world_geo, fill = "grey90")+
  geom_sf(data = pop_sf %>% filter(pop != "MoarBay"), 
          aes(fill = fct_reorder(pop,desc(lat))), color = "black",size = 5,pch = 21)+
  coord_sf(c(-81.75, -76), ylim = c(51.25, 54.75))+
  scale_fill_manual(values = paletteer_d(`"awtools::bpalette"`), name = "Population (n)")+
  scale_size(range=c(0, .5)) +
  annotate(geom = "text", 
           x = -80.7, y = 54, 
           label = "James\nBay", 
           color = "lightblue", 
           size = 8,
           fontface = "italic")+
  ggspatial::annotation_scale(
    location = "tl",
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book"
  ) +
  theme(panel.background = element_rect(fill = "white"),
        axis.text=element_text(size=18),
        axis.title = element_blank(),
        legend.background = element_blank(),
        strip.text.x = element_text(size = 26),
        panel.border = element_rect(colour = "black", fill = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18))


ggsave("~/Desktop/Baie_James_Paper/map_bt.jpeg",width =10, height =10)


