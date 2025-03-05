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
library(maps)
library(mapdata)
library(viridis)
library(ggsn)
library(devtools)
library(ggsn)
library(readxl)


########################
# 2 - Load Data  #
########################

#World map background
world <- map("world", fill=TRUE, plot=FALSE)

world <- ne_countries(
  country = "Canada",
  scale = "large", returnclass = "sf"
)

world_geo <- st_set_crs(world, "EPSG:4326")

#Sampling location data

##Reformatage de long/lat à partir du document original

pop <- read.table("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/pop_sampled.txt",header = T)

pop_bt <- pop %>% filter(!pop_code %in% c("RUP","NOT","JAC"))

pop_sf <- st_as_sf(pop_bt, coords = c("long", "lat"), crs = 4326)

pop_sf$lat <- pop_bt$lat


#-----------#
# James Bay #
#-----------#

#-------------------#
#--- Brook Trout ---#
#-------------------#

rog <- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/roggan/")
kaa <- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/kapsaoui/")
bia <- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/bia/")
lagr <- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/lagrande/")
aqu<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/aquatuc/")
bea<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/beaver/")
sab<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/sabacunica/")
sabl<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/sabacunica_l/")
ren<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/renoyer/")
conn<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/conn/")
eas<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/eastmain/")
jac<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/jack/")
rup<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/rupert/")
bro<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/broadback/")
not<- st_read("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/data/nottaway/")


pop_coords <- as.data.frame(sf::st_coordinates(pop_sf))
pop_coords$pop <- pop_sf$pop

n_pop_bt <- read_xlsx("../metadata.xlsx",sheet = 4)

pop_sf %<>% left_join(n_pop_bt,by = c("pop_code"="site"))

c("#d53e4f",
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

ggplot() +
  geom_sf(data = world_geo,
          fill = "grey85",
          color = "black") +
  geom_sf(data = rog,
          color = "white",
          size = 10) +
  geom_sf(data = kaa, color = "white") +
  geom_sf(data = bia, color = "white") +
  geom_sf(data = lagr, color = "white") +
  geom_sf(data = aqu, color = "white") +
  geom_sf(data = bea, color = "white") +
  geom_sf(data = sab, color = "white") +
  geom_sf(data = sabl, color = "white") +
  geom_sf(data = ren, color = "white") +
  geom_sf(data = conn, color = "white") +
  geom_sf(data = eas, color = "white") +
  geom_sf(data = jac, color = "white") +
  geom_sf(data = rup, color = "white") +
  geom_sf(data = bro, color = "white") +
  geom_sf(data = not, color = "white") +
  geom_sf(
    data = pop_sf %>% filter(!pop %in% c("MoarBay")),
    aes(size = n_wt_miss, fill = reorder(pop_code, lat)),
    color = "black",
    pch = 21
  ) +
  scale_fill_manual(
    values =  c(
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
  ) +
  coord_sf(c(-81.75,-77.5), ylim = c(51.1, 54.75)) +
  ggspatial::annotation_scale(
    location = "tl",
    text_cex = 1,
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book"
  ) +
  north(
    location = "topleft",
    x.min = -81.75,
    x.max = -77,
    y.min = 51.1,
    y.max = 54.75
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 18),
    axis.title = element_blank(),
    strip.text.x = element_text(size = 26),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.background = element_blank(),
    legend.position = c(0.1, 0.75),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  ) +
  guides(fill = "none")


ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/figures/James_Bay_all_sites_bt.png")


#-------------------#
#---- Whitefish ----#
#-------------------#

pop_wf <- pop %>% filter(!pop_code %in% c("BRO","AQU"))

pop_sf_wf <- st_as_sf(pop_wf, coords = c("long", "lat"), crs = 4326)

pop_sf_wf$lat <- pop_wf$lat

n_pop_wf <- read_xlsx("../metadata.xlsx",sheet = 2)

pop_sf_wf %<>% left_join(n_pop_wf,by = c("pop_code"="site"))


ggplot() +
  geom_sf(data = world_geo,
          fill = "grey85",
          color = "black") +
  geom_sf(data = rog,
          color = "white",
          size = 10) +
  geom_sf(data = kaa, color = "white") +
  geom_sf(data = bia, color = "white") +
  geom_sf(data = lagr, color = "white") +
  geom_sf(data = aqu, color = "white") +
  geom_sf(data = bea, color = "white") +
  geom_sf(data = sab, color = "white") +
  geom_sf(data = sabl, color = "white") +
  geom_sf(data = ren, color = "white") +
  geom_sf(data = conn, color = "white") +
  geom_sf(data = eas, color = "white") +
  geom_sf(data = jac, color = "white") +
  geom_sf(data = rup, color = "white") +
  geom_sf(data = bro, color = "white") +
  geom_sf(data = not, color = "white") +
  geom_sf(
    data = pop_sf_wf %>% filter(!pop %in% c("MoarBay")),
    aes(size = n_wt_miss, fill = reorder(pop_code, lat)),
    color = "black",
    pch = 21
  ) +
  scale_fill_manual(values = 
    c("#A50026FF",
      "#DD3D2DFF",
      "#F67E4BFF",
      "#FDB366FF",
      "#FEDA8BFF",
      "#fddbc7",
      "#EAECCCFF",
      "#abdda4",
      "#74add1",
      "#4575b4",
      "#313695",
      "#5e4fa2")
  ) +
  coord_sf(c(-81.75,-77.5), ylim = c(51.1, 54.75)) +
  ggspatial::annotation_scale(
    location = "tl",
    text_cex = 1,
    bar_cols = c("grey60", "white"),
    text_family = "ArcherPro Book"
  ) +
  north(
    location = "topleft",
    x.min = -81.75,
    x.max = -77,
    y.min = 51.1,
    y.max = 54.75
  ) +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 18),
    axis.title = element_blank(),
    strip.text.x = element_text(size = 26),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_text(size = 18),
    legend.text = element_text(size = 14),
    legend.background = element_blank(),
    legend.position = c(0.1, 0.75),
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  ) +
  guides(fill = "none")


ggsave("~/Desktop/Baie_James_Paper/01_population_genomics/comparative_analysis/figures/James_Bay_all_sites_wf.png")

#Whitefish




